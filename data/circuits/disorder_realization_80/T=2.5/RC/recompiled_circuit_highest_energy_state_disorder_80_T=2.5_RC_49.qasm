OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12240527) q[0];
sx q[0];
rz(-0.91949099) q[0];
sx q[0];
rz(-1.1991731) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(-2.5993102) q[1];
sx q[1];
rz(-1.5195001) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8336415) q[0];
sx q[0];
rz(-2.184377) q[0];
sx q[0];
rz(-1.475901) q[0];
rz(2.7678732) q[2];
sx q[2];
rz(-1.9961832) q[2];
sx q[2];
rz(-2.579414) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5120843) q[1];
sx q[1];
rz(-2.0633882) q[1];
sx q[1];
rz(-1.1118481) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24418045) q[3];
sx q[3];
rz(-1.5236519) q[3];
sx q[3];
rz(-2.3650996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2970994) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(-1.7681047) q[2];
rz(0.80396906) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.1871626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74541575) q[0];
sx q[0];
rz(-0.71393037) q[0];
sx q[0];
rz(-0.3748689) q[0];
rz(-0.51775852) q[1];
sx q[1];
rz(-1.9381783) q[1];
sx q[1];
rz(-1.3688603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3387951) q[0];
sx q[0];
rz(-1.5708923) q[0];
sx q[0];
rz(3.1406458) q[0];
rz(2.9451319) q[2];
sx q[2];
rz(-1.5972553) q[2];
sx q[2];
rz(-0.65633869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6102099) q[1];
sx q[1];
rz(-1.0730337) q[1];
sx q[1];
rz(-2.0396292) q[1];
x q[2];
rz(0.13624713) q[3];
sx q[3];
rz(-1.3582503) q[3];
sx q[3];
rz(-3.1177136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.123473) q[2];
sx q[2];
rz(-2.4407083) q[2];
sx q[2];
rz(-2.4269721) q[2];
rz(-2.9361652) q[3];
sx q[3];
rz(-1.3222062) q[3];
sx q[3];
rz(1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6767204) q[0];
sx q[0];
rz(-2.1045852) q[0];
sx q[0];
rz(-2.6439457) q[0];
rz(-0.63703713) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(3.0348437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15539385) q[0];
sx q[0];
rz(-1.2717016) q[0];
sx q[0];
rz(0.22217447) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.006615482) q[2];
sx q[2];
rz(-0.90108904) q[2];
sx q[2];
rz(2.1257328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5046204) q[1];
sx q[1];
rz(-2.2832738) q[1];
sx q[1];
rz(3.0709355) q[1];
rz(-1.2744997) q[3];
sx q[3];
rz(-1.1856286) q[3];
sx q[3];
rz(-0.22814685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60483021) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(-0.19482782) q[2];
rz(3.1365862) q[3];
sx q[3];
rz(-0.61401335) q[3];
sx q[3];
rz(-2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0358589) q[0];
sx q[0];
rz(-2.6073313) q[0];
sx q[0];
rz(1.0400829) q[0];
rz(2.9478574) q[1];
sx q[1];
rz(-1.6400784) q[1];
sx q[1];
rz(0.74660444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0464942) q[0];
sx q[0];
rz(-1.455716) q[0];
sx q[0];
rz(1.1287354) q[0];
x q[1];
rz(0.87073054) q[2];
sx q[2];
rz(-1.680003) q[2];
sx q[2];
rz(2.2094215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0185701) q[1];
sx q[1];
rz(-0.69594068) q[1];
sx q[1];
rz(-3.0607515) q[1];
rz(-pi) q[2];
rz(2.7203619) q[3];
sx q[3];
rz(-2.879749) q[3];
sx q[3];
rz(-0.94294846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99183434) q[2];
sx q[2];
rz(-1.4703625) q[2];
sx q[2];
rz(-2.9370918) q[2];
rz(1.1931984) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(-0.96113718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0123154) q[0];
sx q[0];
rz(-0.17450541) q[0];
sx q[0];
rz(1.0804863) q[0];
rz(2.7642545) q[1];
sx q[1];
rz(-1.5107379) q[1];
sx q[1];
rz(-2.2468755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52496929) q[0];
sx q[0];
rz(-1.7979413) q[0];
sx q[0];
rz(-0.36360111) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65308833) q[2];
sx q[2];
rz(-2.2342367) q[2];
sx q[2];
rz(1.1764248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1354358) q[1];
sx q[1];
rz(-2.4099011) q[1];
sx q[1];
rz(0.77349551) q[1];
rz(0.36578806) q[3];
sx q[3];
rz(-2.180763) q[3];
sx q[3];
rz(-0.47780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6696024) q[2];
sx q[2];
rz(-2.695638) q[2];
sx q[2];
rz(3.0338244) q[2];
rz(-2.9660411) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(0.56041437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046722978) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(1.0119337) q[0];
rz(-1.5049505) q[1];
sx q[1];
rz(-0.71957809) q[1];
sx q[1];
rz(-1.0171657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5115248) q[0];
sx q[0];
rz(-1.8662794) q[0];
sx q[0];
rz(-2.9620671) q[0];
x q[1];
rz(0.58508137) q[2];
sx q[2];
rz(-0.85137109) q[2];
sx q[2];
rz(-2.35308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.519327) q[1];
sx q[1];
rz(-0.19058386) q[1];
sx q[1];
rz(-0.045585338) q[1];
x q[2];
rz(1.8886376) q[3];
sx q[3];
rz(-0.46964619) q[3];
sx q[3];
rz(-1.8152678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4096058) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(2.1448263) q[2];
rz(-3.0766727) q[3];
sx q[3];
rz(-1.7305948) q[3];
sx q[3];
rz(-1.1499278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053059269) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(0.21993318) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(-1.251108) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673888) q[0];
sx q[0];
rz(-2.5944864) q[0];
sx q[0];
rz(-1.0542458) q[0];
x q[1];
rz(1.6701397) q[2];
sx q[2];
rz(-1.2867498) q[2];
sx q[2];
rz(2.5104475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1312853) q[1];
sx q[1];
rz(-0.96988011) q[1];
sx q[1];
rz(2.3405398) q[1];
x q[2];
rz(-0.12507579) q[3];
sx q[3];
rz(-2.2975058) q[3];
sx q[3];
rz(-2.7419326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0327586) q[2];
sx q[2];
rz(-0.83403504) q[2];
sx q[2];
rz(1.533482) q[2];
rz(-1.1533302) q[3];
sx q[3];
rz(-2.0861552) q[3];
sx q[3];
rz(-1.6495033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.15025) q[0];
sx q[0];
rz(-2.7582176) q[0];
sx q[0];
rz(-0.94069329) q[0];
rz(-0.13537814) q[1];
sx q[1];
rz(-1.4043413) q[1];
sx q[1];
rz(-1.0248331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0198284) q[0];
sx q[0];
rz(-1.8786977) q[0];
sx q[0];
rz(2.204049) q[0];
x q[1];
rz(3.054744) q[2];
sx q[2];
rz(-2.5321143) q[2];
sx q[2];
rz(2.0079835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35541818) q[1];
sx q[1];
rz(-0.23705951) q[1];
sx q[1];
rz(2.5393344) q[1];
rz(-0.41127326) q[3];
sx q[3];
rz(-1.6182634) q[3];
sx q[3];
rz(0.91225831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7950874) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-2.4915462) q[2];
rz(0.58642379) q[3];
sx q[3];
rz(-1.4805877) q[3];
sx q[3];
rz(-2.39095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2670249) q[0];
sx q[0];
rz(-2.6331007) q[0];
sx q[0];
rz(-1.9961927) q[0];
rz(1.8264495) q[1];
sx q[1];
rz(-1.2329085) q[1];
sx q[1];
rz(0.68793908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852649) q[0];
sx q[0];
rz(-1.9419799) q[0];
sx q[0];
rz(0.20275499) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97369377) q[2];
sx q[2];
rz(-1.2427689) q[2];
sx q[2];
rz(-2.6231678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8677788) q[1];
sx q[1];
rz(-1.7325337) q[1];
sx q[1];
rz(-2.383197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37738684) q[3];
sx q[3];
rz(-2.0321369) q[3];
sx q[3];
rz(-2.2113706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.410586) q[2];
sx q[2];
rz(-1.1062016) q[2];
sx q[2];
rz(0.86622396) q[2];
rz(2.2335562) q[3];
sx q[3];
rz(-0.76483813) q[3];
sx q[3];
rz(-1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.9662125) q[0];
sx q[0];
rz(-1.9575653) q[0];
sx q[0];
rz(0.41561919) q[0];
rz(0.79187727) q[1];
sx q[1];
rz(-1.917058) q[1];
sx q[1];
rz(2.9143639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450175) q[0];
sx q[0];
rz(-1.554721) q[0];
sx q[0];
rz(-1.5825558) q[0];
x q[1];
rz(-0.80663075) q[2];
sx q[2];
rz(-1.6268332) q[2];
sx q[2];
rz(-0.12412589) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0172248) q[1];
sx q[1];
rz(-1.827888) q[1];
sx q[1];
rz(-2.5460047) q[1];
rz(-pi) q[2];
rz(2.4070508) q[3];
sx q[3];
rz(-0.65207043) q[3];
sx q[3];
rz(1.861524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26433429) q[2];
sx q[2];
rz(-0.72089583) q[2];
sx q[2];
rz(0.43992511) q[2];
rz(0.59709966) q[3];
sx q[3];
rz(-2.4720981) q[3];
sx q[3];
rz(-1.3574903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4933585) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(1.413912) q[1];
sx q[1];
rz(-2.1205817) q[1];
sx q[1];
rz(-3.0243712) q[1];
rz(-1.6261423) q[2];
sx q[2];
rz(-1.9774441) q[2];
sx q[2];
rz(-1.0033506) q[2];
rz(0.24735484) q[3];
sx q[3];
rz(-0.83360278) q[3];
sx q[3];
rz(-3.0616888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
