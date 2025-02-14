OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8799514) q[0];
sx q[0];
rz(-2.9961442) q[0];
sx q[0];
rz(2.8737336) q[0];
rz(1.574006) q[1];
sx q[1];
rz(-2.9730453) q[1];
sx q[1];
rz(-0.57810098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0165839) q[0];
sx q[0];
rz(-1.2947695) q[0];
sx q[0];
rz(-0.87645032) q[0];
x q[1];
rz(2.9801951) q[2];
sx q[2];
rz(-0.70490743) q[2];
sx q[2];
rz(-0.43625956) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0069258) q[1];
sx q[1];
rz(-0.6158661) q[1];
sx q[1];
rz(2.3872972) q[1];
rz(-pi) q[2];
rz(-0.86712305) q[3];
sx q[3];
rz(-1.4438585) q[3];
sx q[3];
rz(2.1849887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1987004) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(-1.4760419) q[2];
rz(-2.5623411) q[3];
sx q[3];
rz(-1.1544635) q[3];
sx q[3];
rz(-2.4756685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2838374) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(-0.055140821) q[0];
rz(2.227123) q[1];
sx q[1];
rz(-1.4417459) q[1];
sx q[1];
rz(-1.2158016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5343842) q[0];
sx q[0];
rz(-1.5522389) q[0];
sx q[0];
rz(0.062848363) q[0];
rz(2.0631172) q[2];
sx q[2];
rz(-1.096306) q[2];
sx q[2];
rz(0.98301393) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7489596) q[1];
sx q[1];
rz(-1.5892423) q[1];
sx q[1];
rz(-0.39876826) q[1];
rz(-pi) q[2];
rz(-1.7082105) q[3];
sx q[3];
rz(-2.0169037) q[3];
sx q[3];
rz(2.844939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2083644) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(-2.5210157) q[2];
rz(1.4536475) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(-2.0645963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2700014) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(-1.510386) q[0];
rz(-0.67391467) q[1];
sx q[1];
rz(-1.2951415) q[1];
sx q[1];
rz(-1.5162226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6166068) q[0];
sx q[0];
rz(-1.5459037) q[0];
sx q[0];
rz(-2.1014433) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7582232) q[2];
sx q[2];
rz(-2.1738833) q[2];
sx q[2];
rz(2.4281339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.080708) q[1];
sx q[1];
rz(-1.9906033) q[1];
sx q[1];
rz(1.8194729) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55530352) q[3];
sx q[3];
rz(-1.8232937) q[3];
sx q[3];
rz(-0.50336526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38500938) q[2];
sx q[2];
rz(-1.5998806) q[2];
sx q[2];
rz(-0.6655244) q[2];
rz(-0.74337983) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(0.22751787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7706364) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(-1.8413405) q[0];
rz(0.11996732) q[1];
sx q[1];
rz(-1.080546) q[1];
sx q[1];
rz(-1.0467451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9366953) q[0];
sx q[0];
rz(-3.129382) q[0];
sx q[0];
rz(-1.3365251) q[0];
rz(1.1566126) q[2];
sx q[2];
rz(-2.3799172) q[2];
sx q[2];
rz(-1.7047395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1018104) q[1];
sx q[1];
rz(-0.48947316) q[1];
sx q[1];
rz(-1.6881949) q[1];
rz(-pi) q[2];
rz(2.1943521) q[3];
sx q[3];
rz(-2.452932) q[3];
sx q[3];
rz(2.6332444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7946502) q[2];
sx q[2];
rz(-1.0618989) q[2];
sx q[2];
rz(-2.3700355) q[2];
rz(-2.0906585) q[3];
sx q[3];
rz(-2.1080878) q[3];
sx q[3];
rz(-1.9482025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13450384) q[0];
sx q[0];
rz(-1.6935885) q[0];
sx q[0];
rz(0.80528468) q[0];
rz(-1.443642) q[1];
sx q[1];
rz(-0.81595683) q[1];
sx q[1];
rz(-0.03646341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591025) q[0];
sx q[0];
rz(-0.56942372) q[0];
sx q[0];
rz(-1.481056) q[0];
x q[1];
rz(-1.6777322) q[2];
sx q[2];
rz(-2.5707939) q[2];
sx q[2];
rz(-0.53174461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.830424) q[1];
sx q[1];
rz(-2.0474381) q[1];
sx q[1];
rz(3.0817506) q[1];
x q[2];
rz(1.7355326) q[3];
sx q[3];
rz(-2.142557) q[3];
sx q[3];
rz(0.71672196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61520758) q[2];
sx q[2];
rz(-1.1247164) q[2];
sx q[2];
rz(2.5649694) q[2];
rz(1.5455101) q[3];
sx q[3];
rz(-0.85448623) q[3];
sx q[3];
rz(-0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312408) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(0.59979576) q[0];
rz(0.80232969) q[1];
sx q[1];
rz(-1.0917412) q[1];
sx q[1];
rz(-0.64782992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14688705) q[0];
sx q[0];
rz(-2.0118666) q[0];
sx q[0];
rz(0.85710454) q[0];
rz(-pi) q[1];
rz(-3.134583) q[2];
sx q[2];
rz(-1.3607549) q[2];
sx q[2];
rz(-1.5633068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5831754) q[1];
sx q[1];
rz(-2.1523471) q[1];
sx q[1];
rz(1.2152332) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69979005) q[3];
sx q[3];
rz(-2.095053) q[3];
sx q[3];
rz(1.4211224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1442147) q[2];
sx q[2];
rz(-0.37313676) q[2];
sx q[2];
rz(-2.0443661) q[2];
rz(1.0002452) q[3];
sx q[3];
rz(-0.92366832) q[3];
sx q[3];
rz(-1.8358811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(-2.9423998) q[0];
rz(-1.8432603) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(-2.4995506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4810774) q[0];
sx q[0];
rz(-0.84781983) q[0];
sx q[0];
rz(-2.8882746) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2175353) q[2];
sx q[2];
rz(-0.63093189) q[2];
sx q[2];
rz(1.0756878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.047028001) q[1];
sx q[1];
rz(-0.21512261) q[1];
sx q[1];
rz(-0.38472979) q[1];
x q[2];
rz(-2.0208259) q[3];
sx q[3];
rz(-1.7535257) q[3];
sx q[3];
rz(-2.8830626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9127427) q[2];
sx q[2];
rz(-2.0759089) q[2];
sx q[2];
rz(1.9672811) q[2];
rz(-2.2187388) q[3];
sx q[3];
rz(-2.703981) q[3];
sx q[3];
rz(0.16930425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6628424) q[0];
sx q[0];
rz(-1.6599382) q[0];
sx q[0];
rz(-2.1424275) q[0];
rz(-2.303458) q[1];
sx q[1];
rz(-0.90845388) q[1];
sx q[1];
rz(2.1160486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2472256) q[0];
sx q[0];
rz(-1.1015838) q[0];
sx q[0];
rz(2.2870025) q[0];
rz(-2.3633358) q[2];
sx q[2];
rz(-2.7822251) q[2];
sx q[2];
rz(0.12675135) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9178114) q[1];
sx q[1];
rz(-1.9635227) q[1];
sx q[1];
rz(-2.1070943) q[1];
x q[2];
rz(2.4225967) q[3];
sx q[3];
rz(-1.9933125) q[3];
sx q[3];
rz(2.7718294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1213558) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(-0.21719246) q[2];
rz(1.1413261) q[3];
sx q[3];
rz(-2.7666028) q[3];
sx q[3];
rz(2.2264437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7104257) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-3.1173832) q[0];
rz(-2.0503893) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(2.3275163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0639125) q[0];
sx q[0];
rz(-1.4850495) q[0];
sx q[0];
rz(-2.4299939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2714839) q[2];
sx q[2];
rz(-0.75557709) q[2];
sx q[2];
rz(2.8062964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7660038) q[1];
sx q[1];
rz(-0.75889041) q[1];
sx q[1];
rz(-1.5689471) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7992448) q[3];
sx q[3];
rz(-1.486393) q[3];
sx q[3];
rz(0.24385246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7476864) q[2];
sx q[2];
rz(-0.86204356) q[2];
sx q[2];
rz(-1.6115335) q[2];
rz(-1.0511506) q[3];
sx q[3];
rz(-1.8007148) q[3];
sx q[3];
rz(3.0401201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53973389) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(2.6644326) q[0];
rz(-1.5513783) q[1];
sx q[1];
rz(-1.662622) q[1];
sx q[1];
rz(1.307055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4068749) q[0];
sx q[0];
rz(-0.85556036) q[0];
sx q[0];
rz(3.0815794) q[0];
rz(-pi) q[1];
x q[1];
rz(2.731852) q[2];
sx q[2];
rz(-1.9671408) q[2];
sx q[2];
rz(0.71045638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8068731) q[1];
sx q[1];
rz(-0.98234896) q[1];
sx q[1];
rz(0.92422331) q[1];
rz(-pi) q[2];
rz(-2.1496592) q[3];
sx q[3];
rz(-2.4483878) q[3];
sx q[3];
rz(2.7424911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9755379) q[2];
sx q[2];
rz(-1.0155948) q[2];
sx q[2];
rz(-2.3892152) q[2];
rz(3.0289529) q[3];
sx q[3];
rz(-2.967716) q[3];
sx q[3];
rz(1.9788474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024121506) q[0];
sx q[0];
rz(-1.385067) q[0];
sx q[0];
rz(-2.0140482) q[0];
rz(2.7680001) q[1];
sx q[1];
rz(-1.5126546) q[1];
sx q[1];
rz(1.0027813) q[1];
rz(2.3898771) q[2];
sx q[2];
rz(-2.5645419) q[2];
sx q[2];
rz(0.6286055) q[2];
rz(1.3704902) q[3];
sx q[3];
rz(-2.3080993) q[3];
sx q[3];
rz(-0.75626683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
