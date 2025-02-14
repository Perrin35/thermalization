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
rz(0.6495629) q[0];
sx q[0];
rz(2.6464638) q[0];
sx q[0];
rz(9.169307) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(4.5402266) q[1];
sx q[1];
rz(11.160688) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509003) q[0];
sx q[0];
rz(-0.003482799) q[0];
sx q[0];
rz(-3.0538959) q[0];
x q[1];
rz(2.6256526) q[2];
sx q[2];
rz(-1.6637515) q[2];
sx q[2];
rz(0.0001212349) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60599323) q[1];
sx q[1];
rz(-0.38315019) q[1];
sx q[1];
rz(-1.5616756) q[1];
rz(-pi) q[2];
rz(-2.5735503) q[3];
sx q[3];
rz(-1.5431817) q[3];
sx q[3];
rz(1.150692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(0.8325038) q[2];
rz(0.80820525) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(2.8830849) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38043624) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(-0.1880745) q[0];
rz(0.07218083) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(-1.6292876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10764194) q[0];
sx q[0];
rz(-2.788262) q[0];
sx q[0];
rz(0.33693846) q[0];
rz(-1.6167655) q[2];
sx q[2];
rz(-1.6013923) q[2];
sx q[2];
rz(1.1083047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.20023274) q[1];
sx q[1];
rz(-1.6299452) q[1];
sx q[1];
rz(0.099353921) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0076804) q[3];
sx q[3];
rz(-0.84656871) q[3];
sx q[3];
rz(-0.40562848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.170681) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(1.2930653) q[2];
rz(-0.34717789) q[3];
sx q[3];
rz(-2.8721589) q[3];
sx q[3];
rz(-0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8186571) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(0.47054189) q[0];
rz(1.200354) q[1];
sx q[1];
rz(-0.73089868) q[1];
sx q[1];
rz(1.2568731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6022105) q[0];
sx q[0];
rz(-0.65927829) q[0];
sx q[0];
rz(-0.97380096) q[0];
x q[1];
rz(-1.7059317) q[2];
sx q[2];
rz(-1.6622512) q[2];
sx q[2];
rz(-2.8460063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4670406) q[1];
sx q[1];
rz(-1.5853154) q[1];
sx q[1];
rz(0.0013703636) q[1];
rz(-0.95945719) q[3];
sx q[3];
rz(-1.3303241) q[3];
sx q[3];
rz(0.38865201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5990126) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-0.92213321) q[2];
rz(1.7982091) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2682997) q[0];
sx q[0];
rz(-2.101185) q[0];
sx q[0];
rz(-2.701395) q[0];
rz(-1.531155) q[1];
sx q[1];
rz(-1.6596158) q[1];
sx q[1];
rz(-0.24756113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.089624) q[0];
sx q[0];
rz(-1.0680833) q[0];
sx q[0];
rz(-0.49552709) q[0];
x q[1];
rz(0.16952408) q[2];
sx q[2];
rz(-1.6598083) q[2];
sx q[2];
rz(1.308987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9369071) q[1];
sx q[1];
rz(-1.8695033) q[1];
sx q[1];
rz(0.077110962) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5400793) q[3];
sx q[3];
rz(-2.2191372) q[3];
sx q[3];
rz(-2.0399486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89746785) q[2];
sx q[2];
rz(-1.030587) q[2];
sx q[2];
rz(-1.6991276) q[2];
rz(2.9407732) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(-2.0012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15631256) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(-0.54541624) q[0];
rz(-3.0975869) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(2.6652179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6355609) q[0];
sx q[0];
rz(-0.2579435) q[0];
sx q[0];
rz(-0.33584612) q[0];
rz(-pi) q[1];
rz(-0.77875359) q[2];
sx q[2];
rz(-2.5980332) q[2];
sx q[2];
rz(-1.5540079) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4551983) q[1];
sx q[1];
rz(-0.87511152) q[1];
sx q[1];
rz(0.2376539) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0073193) q[3];
sx q[3];
rz(-2.0760787) q[3];
sx q[3];
rz(-0.093029417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62006092) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(0.36038348) q[2];
rz(2.9195869) q[3];
sx q[3];
rz(-1.5304151) q[3];
sx q[3];
rz(-1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5510657) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(1.5850413) q[0];
rz(-3.0146154) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(3.051905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2253826) q[0];
sx q[0];
rz(-1.0400335) q[0];
sx q[0];
rz(2.5368774) q[0];
rz(-pi) q[1];
rz(2.3562535) q[2];
sx q[2];
rz(-1.1113104) q[2];
sx q[2];
rz(1.6096514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27881926) q[1];
sx q[1];
rz(-1.3403646) q[1];
sx q[1];
rz(0.50737516) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9305658) q[3];
sx q[3];
rz(-1.0002975) q[3];
sx q[3];
rz(0.24107547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.108532) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(2.8899371) q[2];
rz(-0.039994914) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(0.46372908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43117487) q[0];
sx q[0];
rz(-3.0496821) q[0];
sx q[0];
rz(-2.726626) q[0];
rz(-1.4893432) q[1];
sx q[1];
rz(-0.057561189) q[1];
sx q[1];
rz(0.32589486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1292757) q[0];
sx q[0];
rz(-2.0124448) q[0];
sx q[0];
rz(2.737816) q[0];
x q[1];
rz(1.0524679) q[2];
sx q[2];
rz(-1.4926693) q[2];
sx q[2];
rz(0.6860439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53908097) q[1];
sx q[1];
rz(-1.314114) q[1];
sx q[1];
rz(0.51694444) q[1];
rz(1.8091466) q[3];
sx q[3];
rz(-2.8299667) q[3];
sx q[3];
rz(2.1441604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(1.0464767) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662812) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(0.46324357) q[0];
rz(2.8766368) q[1];
sx q[1];
rz(-0.0015365096) q[1];
sx q[1];
rz(1.6418246) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508598) q[0];
sx q[0];
rz(-1.5895278) q[0];
sx q[0];
rz(2.0035726) q[0];
rz(-pi) q[1];
rz(0.11489377) q[2];
sx q[2];
rz(-0.44821366) q[2];
sx q[2];
rz(-0.047859636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38507515) q[1];
sx q[1];
rz(-0.2236872) q[1];
sx q[1];
rz(-1.7023515) q[1];
rz(-pi) q[2];
rz(-1.8231976) q[3];
sx q[3];
rz(-1.5034672) q[3];
sx q[3];
rz(-2.61781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90193343) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(-1.8233914) q[2];
rz(-0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(1.9759294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5636633) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(0.71459115) q[0];
rz(-0.21358061) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.2108796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7283413) q[0];
sx q[0];
rz(-1.265622) q[0];
sx q[0];
rz(-3.0307709) q[0];
rz(-pi) q[1];
rz(-0.69542747) q[2];
sx q[2];
rz(-1.2992685) q[2];
sx q[2];
rz(0.62715392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80187243) q[1];
sx q[1];
rz(-1.045671) q[1];
sx q[1];
rz(0.32466356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2641422) q[3];
sx q[3];
rz(-0.2444548) q[3];
sx q[3];
rz(2.2742184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96656281) q[2];
sx q[2];
rz(-0.021447072) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(-2.7909732) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(1.9982136) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7570033) q[0];
sx q[0];
rz(-0.92246711) q[0];
sx q[0];
rz(-0.61436999) q[0];
rz(-0.059582926) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(1.7170067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88286582) q[0];
sx q[0];
rz(-2.388991) q[0];
sx q[0];
rz(1.2105788) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5688046) q[2];
sx q[2];
rz(-1.702076) q[2];
sx q[2];
rz(-0.37007331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1150306) q[1];
sx q[1];
rz(-1.6477589) q[1];
sx q[1];
rz(1.413024) q[1];
rz(1.9779102) q[3];
sx q[3];
rz(-1.7686678) q[3];
sx q[3];
rz(-2.1382633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1018579) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-2.6177935) q[2];
rz(1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(-0.4642134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687748) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(-2.4986497) q[2];
sx q[2];
rz(-1.8131154) q[2];
sx q[2];
rz(-2.9396069) q[2];
rz(-0.46661507) q[3];
sx q[3];
rz(-1.3947006) q[3];
sx q[3];
rz(-2.5789921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
