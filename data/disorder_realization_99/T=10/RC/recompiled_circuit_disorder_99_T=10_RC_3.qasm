OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(-0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(2.38148) q[1];
sx q[1];
rz(10.217477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8108609) q[0];
sx q[0];
rz(-2.2415677) q[0];
sx q[0];
rz(-3.1109111) q[0];
x q[1];
rz(0.36926271) q[2];
sx q[2];
rz(-0.57537503) q[2];
sx q[2];
rz(-1.0550261) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3687009) q[1];
sx q[1];
rz(-1.4979616) q[1];
sx q[1];
rz(2.021695) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3000852) q[3];
sx q[3];
rz(-0.91472799) q[3];
sx q[3];
rz(2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(2.7761249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734946) q[0];
sx q[0];
rz(-1.1143648) q[0];
sx q[0];
rz(1.2903851) q[0];
rz(-pi) q[1];
rz(0.89895504) q[2];
sx q[2];
rz(-2.0546753) q[2];
sx q[2];
rz(1.8091786) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47479113) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(2.6067961) q[1];
x q[2];
rz(-0.20110735) q[3];
sx q[3];
rz(-1.6893941) q[3];
sx q[3];
rz(0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6621646) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(-1.1685632) q[0];
x q[1];
rz(0.30655105) q[2];
sx q[2];
rz(-1.2129285) q[2];
sx q[2];
rz(1.4150261) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33774099) q[1];
sx q[1];
rz(-0.9461113) q[1];
sx q[1];
rz(-1.1869663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2080473) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(-2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(1.191656) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062549) q[0];
sx q[0];
rz(-0.65721411) q[0];
sx q[0];
rz(-2.6150136) q[0];
rz(-pi) q[1];
rz(2.4112169) q[2];
sx q[2];
rz(-0.59130284) q[2];
sx q[2];
rz(0.060746047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1034643) q[1];
sx q[1];
rz(-0.80059073) q[1];
sx q[1];
rz(2.994328) q[1];
rz(-pi) q[2];
rz(-1.0563072) q[3];
sx q[3];
rz(-2.0293651) q[3];
sx q[3];
rz(1.6955299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(0.40571037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.273542) q[0];
sx q[0];
rz(-1.2965856) q[0];
sx q[0];
rz(2.4540841) q[0];
rz(-pi) q[1];
rz(2.1359372) q[2];
sx q[2];
rz(-1.7959776) q[2];
sx q[2];
rz(1.856368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5647445) q[1];
sx q[1];
rz(-1.2291359) q[1];
sx q[1];
rz(0.26775743) q[1];
x q[2];
rz(-1.8040854) q[3];
sx q[3];
rz(-2.3217215) q[3];
sx q[3];
rz(1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(-2.3029095) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31081653) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(-2.282417) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(-1.8242594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989527) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(-2.1138666) q[0];
rz(-2.1520734) q[2];
sx q[2];
rz(-0.98810722) q[2];
sx q[2];
rz(-1.6473824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6588905) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(2.3059694) q[1];
rz(-pi) q[2];
rz(2.156593) q[3];
sx q[3];
rz(-1.9337774) q[3];
sx q[3];
rz(-1.7811048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-2.0136925) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.2316661) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8679778) q[0];
sx q[0];
rz(-2.5036739) q[0];
sx q[0];
rz(-2.5173553) q[0];
x q[1];
rz(2.5008051) q[2];
sx q[2];
rz(-1.0377585) q[2];
sx q[2];
rz(-1.2696881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61464804) q[1];
sx q[1];
rz(-1.467418) q[1];
sx q[1];
rz(2.9109216) q[1];
x q[2];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(0.59246078) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97544599) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(-0.59068824) q[0];
rz(-2.5657797) q[2];
sx q[2];
rz(-2.0137219) q[2];
sx q[2];
rz(1.5695614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.2420328) q[1];
sx q[1];
rz(-2.2734103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.153272) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(-1.6944983) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(2.5794199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32414831) q[0];
sx q[0];
rz(-1.9262909) q[0];
sx q[0];
rz(2.4118511) q[0];
rz(-pi) q[1];
x q[1];
rz(0.093900605) q[2];
sx q[2];
rz(-1.0419838) q[2];
sx q[2];
rz(1.4096066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2562981) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(2.2018432) q[1];
x q[2];
rz(0.77769827) q[3];
sx q[3];
rz(-2.5986528) q[3];
sx q[3];
rz(0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-0.047853619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809966) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(3.1050443) q[0];
rz(2.0132952) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(-2.0890582) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1298444) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(-0.5188491) q[1];
x q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3180278) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(0.50280747) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(2.5903493) q[2];
sx q[2];
rz(-1.5521282) q[2];
sx q[2];
rz(2.9350401) q[2];
rz(1.1669284) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];