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
rz(-2.2351216) q[0];
sx q[0];
rz(-1.6989166) q[0];
sx q[0];
rz(2.859681) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(-1.5780916) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1024541) q[0];
sx q[0];
rz(-0.41101563) q[0];
sx q[0];
rz(2.0356112) q[0];
x q[1];
rz(1.5654388) q[2];
sx q[2];
rz(-1.2959576) q[2];
sx q[2];
rz(-1.9071867) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7719745) q[1];
sx q[1];
rz(-1.6965515) q[1];
sx q[1];
rz(2.5096276) q[1];
rz(2.4535577) q[3];
sx q[3];
rz(-1.3322988) q[3];
sx q[3];
rz(-1.9515338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4549183) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(0.61304027) q[2];
rz(-2.6679299) q[3];
sx q[3];
rz(-1.1981755) q[3];
sx q[3];
rz(-1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365725) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(-0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(2.1717333) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365627) q[0];
sx q[0];
rz(-1.9410994) q[0];
sx q[0];
rz(-2.4634393) q[0];
rz(-0.2158176) q[2];
sx q[2];
rz(-1.5310129) q[2];
sx q[2];
rz(1.424274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0326421) q[1];
sx q[1];
rz(-1.3130929) q[1];
sx q[1];
rz(-0.9089246) q[1];
x q[2];
rz(2.0131575) q[3];
sx q[3];
rz(-2.1214888) q[3];
sx q[3];
rz(-0.1503508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2227309) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(-0.53885031) q[2];
rz(3.1239964) q[3];
sx q[3];
rz(-2.9774057) q[3];
sx q[3];
rz(0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(-3.0803296) q[0];
rz(-1.6368658) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(-1.1963074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17249566) q[0];
sx q[0];
rz(-0.75706702) q[0];
sx q[0];
rz(-1.9226546) q[0];
rz(-2.8878651) q[2];
sx q[2];
rz(-2.8710033) q[2];
sx q[2];
rz(1.21278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34431008) q[1];
sx q[1];
rz(-0.6629262) q[1];
sx q[1];
rz(-1.0271038) q[1];
rz(-pi) q[2];
rz(0.34516224) q[3];
sx q[3];
rz(-0.90800873) q[3];
sx q[3];
rz(-1.8927285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4521744) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(2.7090731) q[2];
rz(1.0906667) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(-2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5525621) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(0.20137782) q[0];
rz(0.51678139) q[1];
sx q[1];
rz(-2.7613381) q[1];
sx q[1];
rz(-0.94863272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7761013) q[0];
sx q[0];
rz(-2.1946215) q[0];
sx q[0];
rz(0.43419773) q[0];
rz(-0.8624415) q[2];
sx q[2];
rz(-1.3459599) q[2];
sx q[2];
rz(0.98029691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3624787) q[1];
sx q[1];
rz(-0.50450212) q[1];
sx q[1];
rz(-3.1140559) q[1];
rz(1.473677) q[3];
sx q[3];
rz(-1.0518952) q[3];
sx q[3];
rz(-3.1176381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(-0.55740994) q[2];
rz(-0.080816001) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(-1.935299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4163365) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(-0.15580767) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(-2.9373998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1847398) q[0];
sx q[0];
rz(-0.25190464) q[0];
sx q[0];
rz(-2.8624318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.92856) q[2];
sx q[2];
rz(-1.4430337) q[2];
sx q[2];
rz(3.0054673) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62120092) q[1];
sx q[1];
rz(-2.0546431) q[1];
sx q[1];
rz(-1.3374167) q[1];
rz(-2.3685826) q[3];
sx q[3];
rz(-1.3979619) q[3];
sx q[3];
rz(-2.5484249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2121409) q[2];
sx q[2];
rz(-1.4931623) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(-2.2996969) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(2.1690185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095116422) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(-0.0055775642) q[0];
rz(-1.3439517) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(0.70095789) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537121) q[0];
sx q[0];
rz(-2.0871665) q[0];
sx q[0];
rz(-0.51306458) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0089325) q[2];
sx q[2];
rz(-2.418926) q[2];
sx q[2];
rz(0.36107963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8474104) q[1];
sx q[1];
rz(-0.27833101) q[1];
sx q[1];
rz(-0.14974447) q[1];
x q[2];
rz(-1.0526377) q[3];
sx q[3];
rz(-2.0817882) q[3];
sx q[3];
rz(-2.4571153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8274716) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(1.1673048) q[2];
rz(-2.2589034) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(-2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90671396) q[0];
sx q[0];
rz(-2.0218847) q[0];
sx q[0];
rz(-3.0562905) q[0];
rz(-0.75434297) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(1.0275966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2404382) q[0];
sx q[0];
rz(-1.580997) q[0];
sx q[0];
rz(2.9290694) q[0];
x q[1];
rz(-1.5761887) q[2];
sx q[2];
rz(-0.9741592) q[2];
sx q[2];
rz(2.3062458) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2237146) q[1];
sx q[1];
rz(-0.34668018) q[1];
sx q[1];
rz(-1.9368882) q[1];
rz(-2.6646752) q[3];
sx q[3];
rz(-1.0641885) q[3];
sx q[3];
rz(1.597054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9895642) q[2];
sx q[2];
rz(-2.1542408) q[2];
sx q[2];
rz(2.7222471) q[2];
rz(-2.5069405) q[3];
sx q[3];
rz(-0.80497634) q[3];
sx q[3];
rz(2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.80970508) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(0.69212717) q[0];
rz(2.5634815) q[1];
sx q[1];
rz(-0.41936857) q[1];
sx q[1];
rz(-1.3828166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7958) q[0];
sx q[0];
rz(-1.6383071) q[0];
sx q[0];
rz(-0.019492143) q[0];
x q[1];
rz(-0.67202576) q[2];
sx q[2];
rz(-1.0254854) q[2];
sx q[2];
rz(-1.3915075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18305138) q[1];
sx q[1];
rz(-1.019283) q[1];
sx q[1];
rz(0.27490487) q[1];
rz(2.0973849) q[3];
sx q[3];
rz(-2.3447737) q[3];
sx q[3];
rz(2.6897893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(-1.716506) q[2];
rz(-0.51618451) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(1.9019351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6462964) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(2.7700951) q[0];
rz(-0.83483541) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(3.1239948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3018349) q[0];
sx q[0];
rz(-2.5389414) q[0];
sx q[0];
rz(2.0728803) q[0];
rz(-pi) q[1];
rz(0.73569466) q[2];
sx q[2];
rz(-1.7580877) q[2];
sx q[2];
rz(-0.012030727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.086603348) q[1];
sx q[1];
rz(-1.9412244) q[1];
sx q[1];
rz(2.5309674) q[1];
rz(-pi) q[2];
rz(2.6952088) q[3];
sx q[3];
rz(-2.4842815) q[3];
sx q[3];
rz(3.102234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(0.10031984) q[2];
rz(-2.912168) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(-0.44982287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68596524) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(0.57383865) q[0];
rz(-1.0539184) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(-0.67646772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8343) q[0];
sx q[0];
rz(-2.7678452) q[0];
sx q[0];
rz(0.71889241) q[0];
x q[1];
rz(2.1207379) q[2];
sx q[2];
rz(-0.46560198) q[2];
sx q[2];
rz(2.4684722) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6868848) q[1];
sx q[1];
rz(-1.0446207) q[1];
sx q[1];
rz(1.673102) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0384239) q[3];
sx q[3];
rz(-1.3992953) q[3];
sx q[3];
rz(2.4166783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4544868) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(1.6472316) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(-2.6715265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710707) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(1.8078049) q[1];
sx q[1];
rz(-1.8232657) q[1];
sx q[1];
rz(2.5008536) q[1];
rz(-0.76539466) q[2];
sx q[2];
rz(-1.6366048) q[2];
sx q[2];
rz(-0.66726782) q[2];
rz(-1.3951493) q[3];
sx q[3];
rz(-1.5271389) q[3];
sx q[3];
rz(2.3300119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
