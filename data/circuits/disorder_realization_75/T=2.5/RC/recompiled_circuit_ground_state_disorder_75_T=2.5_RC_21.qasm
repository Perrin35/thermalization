OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.92575443) q[0];
sx q[0];
rz(-0.32045066) q[0];
sx q[0];
rz(3.0132063) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(-0.58712062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7032584) q[0];
sx q[0];
rz(-1.4263784) q[0];
sx q[0];
rz(-0.46972204) q[0];
x q[1];
rz(-0.52337661) q[2];
sx q[2];
rz(-1.9743154) q[2];
sx q[2];
rz(1.425682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93361357) q[1];
sx q[1];
rz(-0.99665239) q[1];
sx q[1];
rz(-2.0962534) q[1];
x q[2];
rz(-1.5774324) q[3];
sx q[3];
rz(-0.38449461) q[3];
sx q[3];
rz(0.16669434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2653653) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(1.5246576) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(2.9558712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8601473) q[0];
sx q[0];
rz(-0.39491072) q[0];
sx q[0];
rz(-1.349378) q[0];
rz(-2.7941864) q[1];
sx q[1];
rz(-1.1888622) q[1];
sx q[1];
rz(-1.4030392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31579298) q[0];
sx q[0];
rz(-0.99195882) q[0];
sx q[0];
rz(1.8301424) q[0];
rz(-0.82282339) q[2];
sx q[2];
rz(-1.1209271) q[2];
sx q[2];
rz(2.8577386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.053005347) q[1];
sx q[1];
rz(-0.8638051) q[1];
sx q[1];
rz(1.7729575) q[1];
rz(-pi) q[2];
rz(0.26936172) q[3];
sx q[3];
rz(-0.90051631) q[3];
sx q[3];
rz(-1.6034077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3226402) q[2];
sx q[2];
rz(-1.8337269) q[2];
sx q[2];
rz(-2.9827706) q[2];
rz(-2.9239376) q[3];
sx q[3];
rz(-1.3889775) q[3];
sx q[3];
rz(-1.7910819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7749216) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(1.2694673) q[0];
rz(2.7216351) q[1];
sx q[1];
rz(-1.0008413) q[1];
sx q[1];
rz(-1.5392083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3184809) q[0];
sx q[0];
rz(-1.2590908) q[0];
sx q[0];
rz(2.0363755) q[0];
rz(-0.69373083) q[2];
sx q[2];
rz(-2.8588534) q[2];
sx q[2];
rz(-3.0764584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2302837) q[1];
sx q[1];
rz(-0.3989987) q[1];
sx q[1];
rz(-1.7085275) q[1];
rz(2.5747803) q[3];
sx q[3];
rz(-1.2897964) q[3];
sx q[3];
rz(-2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0594844) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(-0.19169894) q[2];
rz(-2.5890403) q[3];
sx q[3];
rz(-1.5250165) q[3];
sx q[3];
rz(-1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4100274) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(0.60436526) q[0];
rz(-1.2454237) q[1];
sx q[1];
rz(-2.3912997) q[1];
sx q[1];
rz(2.2818458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5167546) q[0];
sx q[0];
rz(-2.4167912) q[0];
sx q[0];
rz(-1.047931) q[0];
rz(-pi) q[1];
x q[1];
rz(1.672869) q[2];
sx q[2];
rz(-1.5566751) q[2];
sx q[2];
rz(-2.6093687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27450505) q[1];
sx q[1];
rz(-2.9531859) q[1];
sx q[1];
rz(0.26885689) q[1];
x q[2];
rz(2.4212461) q[3];
sx q[3];
rz(-1.7631774) q[3];
sx q[3];
rz(-1.1921574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2802281) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(1.2922618) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(2.0869702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0629405) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(0.31785059) q[0];
rz(1.802313) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(-2.1628765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174216) q[0];
sx q[0];
rz(-1.540543) q[0];
sx q[0];
rz(0.4010648) q[0];
rz(-pi) q[1];
rz(2.9462325) q[2];
sx q[2];
rz(-1.4835412) q[2];
sx q[2];
rz(-1.363014) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.070649466) q[1];
sx q[1];
rz(-1.4448721) q[1];
sx q[1];
rz(0.75448116) q[1];
rz(-2.3581402) q[3];
sx q[3];
rz(-0.26973596) q[3];
sx q[3];
rz(-1.3470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19834441) q[2];
sx q[2];
rz(-0.94503108) q[2];
sx q[2];
rz(-1.1172969) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(-1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4365874) q[0];
sx q[0];
rz(-2.8745108) q[0];
sx q[0];
rz(-1.8319112) q[0];
rz(-2.1197223) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(-1.6530564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40225095) q[0];
sx q[0];
rz(-1.9034804) q[0];
sx q[0];
rz(-2.2646623) q[0];
rz(-pi) q[1];
rz(0.6783046) q[2];
sx q[2];
rz(-1.7540175) q[2];
sx q[2];
rz(1.9187294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1209346) q[1];
sx q[1];
rz(-1.7132732) q[1];
sx q[1];
rz(0.49748904) q[1];
x q[2];
rz(1.009935) q[3];
sx q[3];
rz(-1.6672242) q[3];
sx q[3];
rz(1.4135464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(0.26408163) q[2];
rz(-1.665202) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(-2.7217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61178094) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(2.0773326) q[0];
rz(3.0423959) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(1.681021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997064) q[0];
sx q[0];
rz(-1.0870904) q[0];
sx q[0];
rz(1.4987491) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.869603) q[2];
sx q[2];
rz(-1.5332216) q[2];
sx q[2];
rz(0.36848289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8122892) q[1];
sx q[1];
rz(-1.9871622) q[1];
sx q[1];
rz(-1.8805481) q[1];
rz(-pi) q[2];
rz(1.1292633) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(0.03015524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(0.59949818) q[2];
rz(-0.89764578) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(-3.1035778) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(1.7401485) q[0];
rz(3.0701045) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(-1.00114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18703774) q[0];
sx q[0];
rz(-2.0164549) q[0];
sx q[0];
rz(-1.8513239) q[0];
rz(-pi) q[1];
rz(-0.66364786) q[2];
sx q[2];
rz(-1.2599242) q[2];
sx q[2];
rz(-3.1093017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3197131) q[1];
sx q[1];
rz(-2.1997169) q[1];
sx q[1];
rz(-1.7340388) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97875217) q[3];
sx q[3];
rz(-2.1518937) q[3];
sx q[3];
rz(-1.6421902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1772168) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(-0.56979805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30963323) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(0.69044789) q[0];
rz(-0.18028232) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(2.8188425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.454754) q[0];
sx q[0];
rz(-1.4658064) q[0];
sx q[0];
rz(-2.9872894) q[0];
rz(-pi) q[1];
rz(2.1229136) q[2];
sx q[2];
rz(-1.565287) q[2];
sx q[2];
rz(1.6679279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9195936) q[1];
sx q[1];
rz(-2.1173491) q[1];
sx q[1];
rz(2.2080457) q[1];
rz(0.9154824) q[3];
sx q[3];
rz(-0.99952664) q[3];
sx q[3];
rz(0.20184982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47101578) q[2];
sx q[2];
rz(-1.3479503) q[2];
sx q[2];
rz(-2.1237109) q[2];
rz(-0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(-1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(-0.92593431) q[0];
rz(0.13735859) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(-2.5416809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1178499) q[0];
sx q[0];
rz(-1.568092) q[0];
sx q[0];
rz(0.54393025) q[0];
rz(-0.62905797) q[2];
sx q[2];
rz(-2.0355844) q[2];
sx q[2];
rz(2.128922) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2221189) q[1];
sx q[1];
rz(-1.6012234) q[1];
sx q[1];
rz(-2.1595776) q[1];
x q[2];
rz(2.0447015) q[3];
sx q[3];
rz(-1.2792493) q[3];
sx q[3];
rz(2.4959559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4124734) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(-2.4998383) q[2];
rz(-1.1563835) q[3];
sx q[3];
rz(-2.3803847) q[3];
sx q[3];
rz(1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079035096) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(-1.4398126) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(1.769968) q[2];
sx q[2];
rz(-1.8225852) q[2];
sx q[2];
rz(1.5612623) q[2];
rz(-0.84550459) q[3];
sx q[3];
rz(-0.68544023) q[3];
sx q[3];
rz(-1.1123085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
