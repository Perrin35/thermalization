OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(-2.1332512) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158186) q[0];
sx q[0];
rz(-1.649722) q[0];
sx q[0];
rz(1.1962587) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54097367) q[2];
sx q[2];
rz(-2.7056575) q[2];
sx q[2];
rz(1.774314) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5656149) q[1];
sx q[1];
rz(-1.7276141) q[1];
sx q[1];
rz(0.48139907) q[1];
rz(-pi) q[2];
rz(1.6761682) q[3];
sx q[3];
rz(-0.19427768) q[3];
sx q[3];
rz(-2.7716178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3081554) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(-2.4885139) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.7750875) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.3436226) q[0];
rz(1.9812745) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(2.5210099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337429) q[0];
sx q[0];
rz(-1.7622927) q[0];
sx q[0];
rz(1.2194077) q[0];
rz(-1.1158022) q[2];
sx q[2];
rz(-2.4206851) q[2];
sx q[2];
rz(1.7146595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.008604) q[1];
sx q[1];
rz(-0.78215608) q[1];
sx q[1];
rz(2.6166366) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9602523) q[3];
sx q[3];
rz(-0.6529633) q[3];
sx q[3];
rz(-3.0137526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4429984) q[2];
sx q[2];
rz(-1.5195165) q[2];
sx q[2];
rz(2.911574) q[2];
rz(1.5864141) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(-0.27568278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25552937) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(-1.080876) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(3.0562775) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9832343) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(0.73617022) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5695989) q[2];
sx q[2];
rz(-1.3912918) q[2];
sx q[2];
rz(0.010590503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4115632) q[1];
sx q[1];
rz(-1.7021028) q[1];
sx q[1];
rz(-0.16741411) q[1];
rz(-pi) q[2];
rz(-1.6162698) q[3];
sx q[3];
rz(-0.95004987) q[3];
sx q[3];
rz(2.5705702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(-0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(-1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.008217) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(0.38544449) q[0];
rz(1.7281945) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(2.8916496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88943931) q[0];
sx q[0];
rz(-1.1595386) q[0];
sx q[0];
rz(2.2338339) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5416597) q[2];
sx q[2];
rz(-3.0351808) q[2];
sx q[2];
rz(2.0787042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0562607) q[1];
sx q[1];
rz(-2.0592494) q[1];
sx q[1];
rz(-1.7621098) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79367359) q[3];
sx q[3];
rz(-1.4303615) q[3];
sx q[3];
rz(2.2190203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46127737) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054955) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(-0.48511037) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1331951) q[0];
sx q[0];
rz(-0.48497619) q[0];
sx q[0];
rz(2.4147245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9747462) q[2];
sx q[2];
rz(-1.8649668) q[2];
sx q[2];
rz(2.6427513) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5768535) q[1];
sx q[1];
rz(-0.67138636) q[1];
sx q[1];
rz(-1.6480084) q[1];
rz(-pi) q[2];
rz(-2.740871) q[3];
sx q[3];
rz(-1.6288789) q[3];
sx q[3];
rz(-2.7207014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(-2.2122993) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-2.3727681) q[0];
rz(1.7898412) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(2.2519055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79737299) q[0];
sx q[0];
rz(-1.780172) q[0];
sx q[0];
rz(1.1816417) q[0];
rz(-1.5264411) q[2];
sx q[2];
rz(-1.5694478) q[2];
sx q[2];
rz(1.9951472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5714214) q[1];
sx q[1];
rz(-1.4246877) q[1];
sx q[1];
rz(3.094638) q[1];
x q[2];
rz(-1.4558639) q[3];
sx q[3];
rz(-2.0191666) q[3];
sx q[3];
rz(-1.861972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49030226) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5550522) q[2];
rz(-1.8709315) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(0.13203013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.1726058) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(-0.83089685) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(2.3544618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68987304) q[0];
sx q[0];
rz(-1.4362808) q[0];
sx q[0];
rz(-3.0089507) q[0];
rz(2.2550341) q[2];
sx q[2];
rz(-2.0066924) q[2];
sx q[2];
rz(-0.8800216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5097376) q[1];
sx q[1];
rz(-0.69586772) q[1];
sx q[1];
rz(-0.79001714) q[1];
rz(-pi) q[2];
rz(-0.55946405) q[3];
sx q[3];
rz(-1.9064404) q[3];
sx q[3];
rz(1.0192724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1181011) q[2];
sx q[2];
rz(-0.96994895) q[2];
sx q[2];
rz(3.1084295) q[2];
rz(1.5983332) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(1.6395114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.24932662) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(-1.7484885) q[0];
rz(2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(-1.1005864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7129214) q[0];
sx q[0];
rz(-1.5690737) q[0];
sx q[0];
rz(-0.07600204) q[0];
x q[1];
rz(0.69533657) q[2];
sx q[2];
rz(-2.0954663) q[2];
sx q[2];
rz(-2.6161043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3708027) q[1];
sx q[1];
rz(-2.2375771) q[1];
sx q[1];
rz(-0.043170269) q[1];
rz(1.4310775) q[3];
sx q[3];
rz(-1.9587687) q[3];
sx q[3];
rz(-3.1041077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.042772375) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(2.8501046) q[0];
rz(-1.905722) q[1];
sx q[1];
rz(-1.196685) q[1];
sx q[1];
rz(3.018697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6345469) q[0];
sx q[0];
rz(-1.6222008) q[0];
sx q[0];
rz(2.794751) q[0];
rz(-pi) q[1];
rz(2.9192186) q[2];
sx q[2];
rz(-1.8463328) q[2];
sx q[2];
rz(1.0433407) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0258515) q[1];
sx q[1];
rz(-2.9065955) q[1];
sx q[1];
rz(-2.0629289) q[1];
rz(-1.5294262) q[3];
sx q[3];
rz(-0.14474104) q[3];
sx q[3];
rz(-0.54171692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5128532) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(-0.17744803) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(-1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43750957) q[0];
sx q[0];
rz(-0.45314416) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-2.5240555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3445471) q[0];
sx q[0];
rz(-2.8918242) q[0];
sx q[0];
rz(-2.0956371) q[0];
rz(1.0087291) q[2];
sx q[2];
rz(-0.8265087) q[2];
sx q[2];
rz(1.1240608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3869218) q[1];
sx q[1];
rz(-2.6082509) q[1];
sx q[1];
rz(1.6120595) q[1];
rz(-pi) q[2];
rz(-2.6838949) q[3];
sx q[3];
rz(-0.62333306) q[3];
sx q[3];
rz(1.8231572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2403229) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(2.0591002) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33901535) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-0.045724178) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(0.5068266) q[2];
sx q[2];
rz(-2.0438556) q[2];
sx q[2];
rz(0.023509916) q[2];
rz(2.4847538) q[3];
sx q[3];
rz(-1.4359063) q[3];
sx q[3];
rz(-2.0834854) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
