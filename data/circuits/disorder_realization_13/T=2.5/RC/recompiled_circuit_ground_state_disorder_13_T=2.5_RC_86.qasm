OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(2.9958041) q[0];
sx q[0];
rz(11.267405) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(4.4823449) q[1];
sx q[1];
rz(9.525099) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.318905) q[0];
sx q[0];
rz(-1.5189369) q[0];
sx q[0];
rz(0.6619617) q[0];
rz(-pi) q[1];
rz(-2.6703628) q[2];
sx q[2];
rz(-1.8719851) q[2];
sx q[2];
rz(0.1567086) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24619094) q[1];
sx q[1];
rz(-2.7936087) q[1];
sx q[1];
rz(-0.65234185) q[1];
rz(-pi) q[2];
rz(2.6907608) q[3];
sx q[3];
rz(-2.5019162) q[3];
sx q[3];
rz(1.0420639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7413062) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-2.3347704) q[2];
rz(-2.412879) q[3];
sx q[3];
rz(-0.7586793) q[3];
sx q[3];
rz(-1.2046643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5173986) q[0];
sx q[0];
rz(-2.875858) q[0];
sx q[0];
rz(2.8345795) q[0];
rz(-1.864805) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(1.101864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4808896) q[0];
sx q[0];
rz(-2.1212949) q[0];
sx q[0];
rz(-0.20396978) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5069783) q[2];
sx q[2];
rz(-1.391727) q[2];
sx q[2];
rz(-3.0916391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9872322) q[1];
sx q[1];
rz(-1.4353961) q[1];
sx q[1];
rz(-1.3920648) q[1];
rz(2.3526221) q[3];
sx q[3];
rz(-2.1405102) q[3];
sx q[3];
rz(2.8017442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.028194204) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(-0.096435189) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518799) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(2.7413947) q[0];
rz(-1.6290889) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(2.8834744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982957) q[0];
sx q[0];
rz(-0.27320293) q[0];
sx q[0];
rz(2.1092723) q[0];
rz(-pi) q[1];
rz(3.0417241) q[2];
sx q[2];
rz(-1.6503008) q[2];
sx q[2];
rz(-0.85209633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8761354) q[1];
sx q[1];
rz(-1.5578414) q[1];
sx q[1];
rz(1.5658246) q[1];
rz(-pi) q[2];
rz(2.9003224) q[3];
sx q[3];
rz(-1.4033969) q[3];
sx q[3];
rz(2.3679855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(-0.26767996) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(-1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716008) q[0];
sx q[0];
rz(-1.6331693) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(-3.1056504) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(2.8004004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55732176) q[0];
sx q[0];
rz(-2.1958618) q[0];
sx q[0];
rz(-1.2465501) q[0];
x q[1];
rz(1.8282537) q[2];
sx q[2];
rz(-0.99960589) q[2];
sx q[2];
rz(-1.5549623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2349639) q[1];
sx q[1];
rz(-2.2412657) q[1];
sx q[1];
rz(-2.5690298) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5248923) q[3];
sx q[3];
rz(-0.69717191) q[3];
sx q[3];
rz(-2.8215849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7146032) q[2];
sx q[2];
rz(-2.0689071) q[2];
sx q[2];
rz(-1.535447) q[2];
rz(-2.8793907) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(-0.95482993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75761211) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(2.0641932) q[0];
rz(2.7491838) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(-2.1108625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8547671) q[0];
sx q[0];
rz(-2.1041901) q[0];
sx q[0];
rz(1.0817762) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2472146) q[2];
sx q[2];
rz(-2.6204118) q[2];
sx q[2];
rz(-0.22400907) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2490152) q[1];
sx q[1];
rz(-1.7356145) q[1];
sx q[1];
rz(-1.7958162) q[1];
x q[2];
rz(-1.0038336) q[3];
sx q[3];
rz(-0.58708411) q[3];
sx q[3];
rz(-2.7239012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82421676) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(-2.2281036) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(2.998108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(1.6827253) q[0];
sx q[0];
rz(-2.904628) q[0];
sx q[0];
rz(1.4759901) q[0];
rz(-2.7492375) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(-0.59142339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3680843) q[0];
sx q[0];
rz(-0.71660935) q[0];
sx q[0];
rz(0.98978292) q[0];
rz(-0.5742091) q[2];
sx q[2];
rz(-0.32445766) q[2];
sx q[2];
rz(2.7221219) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6900031) q[1];
sx q[1];
rz(-2.0760447) q[1];
sx q[1];
rz(0.075376781) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2675681) q[3];
sx q[3];
rz(-1.8667779) q[3];
sx q[3];
rz(-0.10914739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(-0.57233125) q[2];
rz(0.22219292) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-0.40827665) q[0];
rz(-0.73221842) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(2.8439723) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0746411) q[0];
sx q[0];
rz(-2.5805001) q[0];
sx q[0];
rz(2.448161) q[0];
x q[1];
rz(-1.5333129) q[2];
sx q[2];
rz(-0.49502326) q[2];
sx q[2];
rz(-1.9062454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.630328) q[1];
sx q[1];
rz(-1.9115834) q[1];
sx q[1];
rz(-2.9326669) q[1];
x q[2];
rz(0.70370378) q[3];
sx q[3];
rz(-1.6856442) q[3];
sx q[3];
rz(-0.59806693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1709661) q[2];
sx q[2];
rz(-1.8793224) q[2];
sx q[2];
rz(-0.69620281) q[2];
rz(-0.89037406) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(-1.4674998) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-0.21275511) q[0];
rz(2.6720324) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(-2.3874217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045005) q[0];
sx q[0];
rz(-1.8593328) q[0];
sx q[0];
rz(0.19591422) q[0];
rz(-pi) q[1];
rz(2.5905072) q[2];
sx q[2];
rz(-0.37691016) q[2];
sx q[2];
rz(-2.2152679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9355802) q[1];
sx q[1];
rz(-1.6846865) q[1];
sx q[1];
rz(-0.64958944) q[1];
x q[2];
rz(-2.2175104) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(-1.9677066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.135005) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(2.3593486) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924454) q[0];
sx q[0];
rz(-0.47098422) q[0];
sx q[0];
rz(-2.1771722) q[0];
rz(-1.2760705) q[1];
sx q[1];
rz(-1.7274906) q[1];
sx q[1];
rz(1.5020348) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3827189) q[0];
sx q[0];
rz(-1.3387464) q[0];
sx q[0];
rz(2.9358572) q[0];
x q[1];
rz(1.1310857) q[2];
sx q[2];
rz(-2.2520116) q[2];
sx q[2];
rz(1.1352254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.561475) q[1];
sx q[1];
rz(-2.0881542) q[1];
sx q[1];
rz(-3.1151732) q[1];
x q[2];
rz(0.067599452) q[3];
sx q[3];
rz(-1.7177594) q[3];
sx q[3];
rz(0.33892469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2532578) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(-2.1962732) q[2];
rz(2.3156598) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(-0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.076040529) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(-2.3384576) q[0];
rz(1.5638634) q[1];
sx q[1];
rz(-1.6604661) q[1];
sx q[1];
rz(-0.28958431) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11293791) q[0];
sx q[0];
rz(-1.5585596) q[0];
sx q[0];
rz(-1.4652246) q[0];
rz(-pi) q[1];
rz(-2.3371465) q[2];
sx q[2];
rz(-2.7000393) q[2];
sx q[2];
rz(-1.785977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3037422) q[1];
sx q[1];
rz(-1.2134064) q[1];
sx q[1];
rz(2.9010495) q[1];
rz(2.8554162) q[3];
sx q[3];
rz(-1.8916142) q[3];
sx q[3];
rz(2.8367354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3760066) q[2];
sx q[2];
rz(-0.12066081) q[2];
sx q[2];
rz(2.181459) q[2];
rz(-2.5586186) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.690602) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(3.1008537) q[1];
sx q[1];
rz(-2.4650885) q[1];
sx q[1];
rz(-3.0104641) q[1];
rz(-2.0410983) q[2];
sx q[2];
rz(-1.294877) q[2];
sx q[2];
rz(-0.69653947) q[2];
rz(-3.0906755) q[3];
sx q[3];
rz(-1.6418727) q[3];
sx q[3];
rz(1.4486936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
