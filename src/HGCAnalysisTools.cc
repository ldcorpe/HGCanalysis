#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include <vector>

using namespace std;

//
float getLambdaForHGCLayer(int hit_layer)
{
  if (hit_layer==1)                       return 0.01; 
  else if(hit_layer>=2 && hit_layer<=11)  return 0.036;
  else if(hit_layer>=12 && hit_layer<=21) return 0.043;
  else if(hit_layer>=22 && hit_layer<=30) return 0.056;
  else if(hit_layer==31)                  return 0.338;
  else if(hit_layer>=32 && hit_layer<=42) return 0.273;
  else if(hit_layer>42)                   return 0.475;
  return 0;
}
//
G4InteractionPositionInfo getInteractionPositionLC(const std::vector<SimTrack> *SimTk, 
						 const std::vector<SimVertex> *SimVtx, float pt )
						//, reco::GenParticle *gen)
{
  // Declare the variables we will need
	G4InteractionPositionInfo toRet; // the reference of the photon conversion vertex, if there is one.
	toRet.pos=math::XYZVectorD(0,0,0); // position of the conversion. if none, return 0.
	toRet.info=0; 
	int trackID = -1; // track ID of genparticle.
	int rawTkMult(0),eTkMult(0); //track mutliplicities, specifically electrons.
	
	// first loop over the available tracks and find the track mathcing the genPhoton. We match using pt, which should be IDENTICAL ebtween gen particle and simtrack.
	for (const SimTrack &vtxTk : *SimTk)
	{

		int tkType=vtxTk.type(); // gives pdg id of track.

		if (abs(tkType)!= 22) continue; // only want photons.
		if( abs(pt - vtxTk.momentum().pt()) > 0.1) continue; // specifically, photons which have the same pt as the genphoton.

	//	std::cout << "[debug ] track " << rawTkMult << ", type " << tkType << ", tk pt " << vtxTk.momentum().pt() << ", true pt " << pt << std::endl;
		trackID = vtxTk.trackId() ; // save the ID of the track corresponding to the gen particle.
		break; // found it? break teh loop.
	}

	//loop over vertices  to find vertex which comes from teh track we found above.
	for (const SimVertex &simVtx : *SimVtx) 
	{

		//require the parent to be the track we found above
		int pIdx( simVtx.parentIndex() );
		if( pIdx!=trackID) continue;  // trackID acts as teh barcode.

		//std::cout << "[debug ] found vertex at x" << simVtx.position().X() << ", y " << simVtx.position().Y() <<", z " << simVtx.position().Z()<< std::endl;

		int vtxIdx(simVtx.vertexId());// ..save the ID of the vertex which has the track as parent.
		
		// now loop over tracks coming from this vertex: look for electrons
		for (const SimTrack &vtxTk : *SimTk)
		{
			int tkVtxIdx( vtxTk.vertIndex() ); 
			if(tkVtxIdx!=vtxIdx) continue; // only look at tracks coming from the above vertex, which is associated with teh genParticle track.

			int tkType=vtxTk.type(); // get PDG id

			//	if (abs(tkType)!= 11) continue;

			rawTkMult++;
				eTkMult      += (abs(tkType)==11); //count number of electrons.
			//std::cout << "[debug ] track " << rawTkMult << ", type " << tkType <<", tk vtx id " << tkVtxIdx<< std::endl;
			//break;
		}

		if(rawTkMult<2) continue; //should have at least one track! 

		toRet.pos=math::XYZVectorD(simVtx.position()); // get position
	//	break;
//		if(rawTkMult==3 && nTkMult==1 && nucleiTkMult==1 && pTkMult==1) toRet.info=1; // not sure what this does.
		}

		return toRet;
	}

	//
	G4InteractionPositionInfo getInteractionPosition(const std::vector<SimTrack> *SimTk, 
			const std::vector<SimVertex> *SimVtx, 
			int barcode)
	{
		G4InteractionPositionInfo toRet;
		toRet.pos=math::XYZVectorD(0,0,0);
		toRet.info=0;

		//loop over vertices
		for (const SimVertex &simVtx : *SimVtx) 
		{
			//require the parent to be the given barcode
			bool noParent( simVtx.noParent() );
			if(noParent) continue;
			int pIdx( simVtx.parentIndex() );
			if( pIdx!=barcode) continue;

			std::cout << "[debug ] found vertex at x" << simVtx.position().X() << ", y " << simVtx.position().Y() <<", z " << simVtx.position().Z()<< std::endl;

			int vtxIdx(simVtx.vertexId());
			int rawTkMult(0),eTkMult(0),gTkMult(0),nTkMult(0),nucleiTkMult(0),pTkMult(0);
			for (const SimTrack &vtxTk : *SimTk)
			{
				int tkVtxIdx( vtxTk.vertIndex() ); 
				if(tkVtxIdx!=vtxIdx) continue;

				int tkType=vtxTk.type();

				//	if (abs(tkType)!= 11) continue;

				rawTkMult++;
				eTkMult      += (abs(tkType)==11);
				gTkMult      += (abs(tkType)==22);
				nTkMult      += (abs(tkType)==2112 || abs(tkType)==2212);
				nucleiTkMult += (abs(tkType)>1000000000);
				pTkMult      += (abs(tkType)==211 || abs(tkType)==111);
				std::cout << "[debug ] track " << rawTkMult << ", type " << tkType <<", tk vtx id " << tkVtxIdx<< std::endl;
			}

			if(rawTkMult<2) continue;

			toRet.pos=math::XYZVectorD(simVtx.position());
			if(rawTkMult==3 && nTkMult==1 && nucleiTkMult==1 && pTkMult==1) toRet.info=1;
			return toRet;
		}

		return toRet;
	}


	//
	std::pair<float,float> getEffSigma(RooRealVar *var, RooAbsPdf *pdf, float wmin,float wmax, float step, float epsilon)
	{
		//get cdf points
		RooAbsReal *cdf = pdf->createCdf(RooArgList(*var));
		float point=wmin;
		vector<pair<float,float> > points;
		while (point <= wmax){
			var->setVal(point);
			if (pdf->getVal() > epsilon){
				points.push_back(pair<float,float>(point,cdf->getVal()));
			}
			point+=step;
		}

		float low = wmin;
		float high = wmax;
		float width = wmax-wmin;
		for (unsigned int i=0; i<points.size(); i++){
			for (unsigned int j=i; j<points.size(); j++){
				float wy = points[j].second - points[i].second;
				if (TMath::Abs(wy-0.683) < epsilon){
					float wx = points[j].first - points[i].first;
					if (wx < width){
						low = points[i].first;
						high = points[j].first;
						width=wx;
					}
				}
			}
		}
		pair<float,float> result(low,high);
		return result;
	}
