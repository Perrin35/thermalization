OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7497082) q[0];
sx q[0];
rz(-2.9449129) q[0];
sx q[0];
rz(-1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084378622) q[0];
sx q[0];
rz(-2.9998261) q[0];
sx q[0];
rz(1.0300107) q[0];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(-1.4717799) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37669668) q[1];
sx q[1];
rz(-0.88250676) q[1];
sx q[1];
rz(0.98801686) q[1];
x q[2];
rz(1.7114867) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(-1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-2.8536318) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156793) q[0];
sx q[0];
rz(-0.26198146) q[0];
sx q[0];
rz(-2.322305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9777074) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(1.5175982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60625633) q[1];
sx q[1];
rz(-1.9793946) q[1];
sx q[1];
rz(2.4819083) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6091299) q[3];
sx q[3];
rz(-1.9209071) q[3];
sx q[3];
rz(2.4083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(1.7193517) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-2.1898988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6355977) q[0];
sx q[0];
rz(-1.6615168) q[0];
sx q[0];
rz(1.3784301) q[0];
rz(2.0864262) q[2];
sx q[2];
rz(-1.7679169) q[2];
sx q[2];
rz(0.65161639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.051085) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53593105) q[3];
sx q[3];
rz(-1.0169573) q[3];
sx q[3];
rz(-1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13016985) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0043871) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(2.7688162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16762776) q[0];
sx q[0];
rz(-2.3229257) q[0];
sx q[0];
rz(1.7276093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77809019) q[2];
sx q[2];
rz(-1.2568682) q[2];
sx q[2];
rz(2.2997466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2671632) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(2.8664385) q[1];
x q[2];
rz(-2.0471441) q[3];
sx q[3];
rz(-1.8653449) q[3];
sx q[3];
rz(-1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(-2.2367031) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596817) q[0];
sx q[0];
rz(-2.7026183) q[0];
sx q[0];
rz(2.3460593) q[0];
rz(0.34110951) q[2];
sx q[2];
rz(-2.038539) q[2];
sx q[2];
rz(1.9175921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91671645) q[1];
sx q[1];
rz(-1.6460878) q[1];
sx q[1];
rz(-3.0776943) q[1];
rz(-2.9910473) q[3];
sx q[3];
rz(-2.4171722) q[3];
sx q[3];
rz(2.5151593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-2.9980998) q[2];
rz(-1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(-1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(1.3751078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28778827) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(0.35065035) q[0];
rz(-pi) q[1];
rz(2.3926922) q[2];
sx q[2];
rz(-1.6237215) q[2];
sx q[2];
rz(-0.63268328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.314234) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(0.41008653) q[1];
x q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-2.0241963) q[3];
sx q[3];
rz(-2.6231433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(-0.14870816) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-2.9490411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(0.08392863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.98691) q[0];
sx q[0];
rz(-1.0202408) q[0];
sx q[0];
rz(2.8793392) q[0];
rz(2.7883045) q[2];
sx q[2];
rz(-1.366426) q[2];
sx q[2];
rz(2.3452961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0107207) q[1];
sx q[1];
rz(-2.2166703) q[1];
sx q[1];
rz(-2.3400397) q[1];
rz(-0.84079068) q[3];
sx q[3];
rz(-1.519324) q[3];
sx q[3];
rz(1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(0.54523462) q[2];
rz(-0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(-0.94582311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686533) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(-0.22139876) q[0];
x q[1];
rz(0.85407599) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(-2.9727109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3555803) q[1];
sx q[1];
rz(-1.5994206) q[1];
sx q[1];
rz(2.7629258) q[1];
rz(0.010057851) q[3];
sx q[3];
rz(-0.66614449) q[3];
sx q[3];
rz(0.87107623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(2.5659134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3424123) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(-0.70223017) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11328463) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(-1.726113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8262144) q[1];
sx q[1];
rz(-1.4634906) q[1];
sx q[1];
rz(-0.025749287) q[1];
rz(-pi) q[2];
rz(0.18746312) q[3];
sx q[3];
rz(-1.0378569) q[3];
sx q[3];
rz(0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(0.021961948) q[2];
rz(0.17045505) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(0.6341933) q[0];
rz(-0.028907396) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-2.172487) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076981912) q[0];
sx q[0];
rz(-1.6425898) q[0];
sx q[0];
rz(-3.0781156) q[0];
rz(-pi) q[1];
rz(-2.3949941) q[2];
sx q[2];
rz(-1.714163) q[2];
sx q[2];
rz(-1.6105086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58418729) q[1];
sx q[1];
rz(-2.4483042) q[1];
sx q[1];
rz(0.28179817) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5892331) q[3];
sx q[3];
rz(-0.5061572) q[3];
sx q[3];
rz(1.6278704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-0.36007145) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(-2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(0.85514851) q[2];
sx q[2];
rz(-2.6519041) q[2];
sx q[2];
rz(2.3408163) q[2];
rz(1.2561856) q[3];
sx q[3];
rz(-1.1369858) q[3];
sx q[3];
rz(2.9997957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
