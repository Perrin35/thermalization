OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3992231) q[0];
sx q[0];
rz(-2.7288781) q[0];
sx q[0];
rz(-0.8362008) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.598624) q[0];
sx q[0];
rz(-2.7592139) q[0];
sx q[0];
rz(-1.3578869) q[0];
x q[1];
rz(-1.8062334) q[2];
sx q[2];
rz(-1.2004235) q[2];
sx q[2];
rz(0.78200227) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5656149) q[1];
sx q[1];
rz(-1.7276141) q[1];
sx q[1];
rz(-2.6601936) q[1];
x q[2];
rz(1.6761682) q[3];
sx q[3];
rz(-0.19427768) q[3];
sx q[3];
rz(-2.7716178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8334373) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(1.0623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7750875) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(1.9812745) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(2.5210099) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68408332) q[0];
sx q[0];
rz(-2.7433222) q[0];
sx q[0];
rz(2.0837559) q[0];
x q[1];
rz(0.90262516) q[2];
sx q[2];
rz(-1.2764954) q[2];
sx q[2];
rz(0.49605344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.008604) q[1];
sx q[1];
rz(-0.78215608) q[1];
sx q[1];
rz(0.52495606) q[1];
x q[2];
rz(-0.28264265) q[3];
sx q[3];
rz(-2.1676554) q[3];
sx q[3];
rz(0.34908906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6985942) q[2];
sx q[2];
rz(-1.5195165) q[2];
sx q[2];
rz(-2.911574) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(-2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25552937) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(1.080876) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(-0.08531514) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9832343) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(2.4054224) q[0];
rz(2.8181452) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(-1.310854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7300295) q[1];
sx q[1];
rz(-1.7021028) q[1];
sx q[1];
rz(-0.16741411) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62123589) q[3];
sx q[3];
rz(-1.6077822) q[3];
sx q[3];
rz(1.0262352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9366511) q[2];
sx q[2];
rz(-0.0095657883) q[2];
sx q[2];
rz(-2.9452475) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(2.096368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008217) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(0.38544449) q[0];
rz(-1.7281945) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(-0.24994303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7628544) q[0];
sx q[0];
rz(-2.1702499) q[0];
sx q[0];
rz(-2.6361639) q[0];
rz(-pi) q[1];
rz(3.0536575) q[2];
sx q[2];
rz(-1.5107949) q[2];
sx q[2];
rz(-2.0363925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6648207) q[1];
sx q[1];
rz(-0.52174924) q[1];
sx q[1];
rz(0.34362766) q[1];
rz(-pi) q[2];
rz(-0.79367359) q[3];
sx q[3];
rz(-1.7112312) q[3];
sx q[3];
rz(-0.92257231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46127737) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(0.2363905) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(-2.4642956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054955) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(2.6564823) q[0];
rz(1.1803892) q[1];
sx q[1];
rz(-1.5943269) q[1];
sx q[1];
rz(2.3050883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0083975365) q[0];
sx q[0];
rz(-0.48497619) q[0];
sx q[0];
rz(2.4147245) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91412394) q[2];
sx q[2];
rz(-0.49490041) q[2];
sx q[2];
rz(2.6656541) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0750351) q[1];
sx q[1];
rz(-1.5227942) q[1];
sx q[1];
rz(2.2407299) q[1];
x q[2];
rz(2.740871) q[3];
sx q[3];
rz(-1.6288789) q[3];
sx q[3];
rz(2.7207014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(-2.3720429) q[2];
rz(2.2122993) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-2.3727681) q[0];
rz(1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(-2.2519055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8990435) q[0];
sx q[0];
rz(-2.7022323) q[0];
sx q[0];
rz(-1.0602632) q[0];
rz(-pi) q[1];
rz(-1.6011997) q[2];
sx q[2];
rz(-3.0972169) q[2];
sx q[2];
rz(0.45472431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8823679) q[1];
sx q[1];
rz(-2.9881757) q[1];
sx q[1];
rz(1.8795561) q[1];
x q[2];
rz(-1.6857288) q[3];
sx q[3];
rz(-1.1224261) q[3];
sx q[3];
rz(1.2796206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49030226) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(-1.8709315) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(0.13203013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1726058) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(1.9975115) q[0];
rz(-2.3106958) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(0.78713083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09296552) q[0];
sx q[0];
rz(-2.9529609) q[0];
sx q[0];
rz(-2.3446871) q[0];
x q[1];
rz(-2.2550341) q[2];
sx q[2];
rz(-1.1349003) q[2];
sx q[2];
rz(2.2615711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28879189) q[1];
sx q[1];
rz(-2.0388985) q[1];
sx q[1];
rz(-2.1062984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5601117) q[3];
sx q[3];
rz(-2.498501) q[3];
sx q[3];
rz(-3.0742925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-3.1084295) q[2];
rz(-1.5983332) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.892266) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(-1.3931042) q[0];
rz(-0.45685592) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(-2.0410062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2835866) q[0];
sx q[0];
rz(-1.6467983) q[0];
sx q[0];
rz(-1.572524) q[0];
rz(-2.216711) q[2];
sx q[2];
rz(-0.98306984) q[2];
sx q[2];
rz(0.6494259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9683084) q[1];
sx q[1];
rz(-1.6047162) q[1];
sx q[1];
rz(0.90356253) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7105152) q[3];
sx q[3];
rz(-1.1828239) q[3];
sx q[3];
rz(0.037484976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0988203) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(-1.5607321) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(-2.1055351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7324657) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(1.2358707) q[1];
sx q[1];
rz(-1.196685) q[1];
sx q[1];
rz(3.018697) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0964121) q[0];
sx q[0];
rz(-1.9171606) q[0];
sx q[0];
rz(-1.5161432) q[0];
rz(-pi) q[1];
rz(0.22237402) q[2];
sx q[2];
rz(-1.2952598) q[2];
sx q[2];
rz(-2.0982519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.025561573) q[1];
sx q[1];
rz(-1.4605547) q[1];
sx q[1];
rz(1.7787537) q[1];
x q[2];
rz(1.5294262) q[3];
sx q[3];
rz(-2.9968516) q[3];
sx q[3];
rz(2.5998757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62873944) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(2.9641446) q[3];
sx q[3];
rz(-0.95824233) q[3];
sx q[3];
rz(-1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040831) q[0];
sx q[0];
rz(-0.45314416) q[0];
sx q[0];
rz(1.7105239) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(2.5240555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7970456) q[0];
sx q[0];
rz(-0.24976845) q[0];
sx q[0];
rz(1.0459556) q[0];
rz(-pi) q[1];
rz(2.3138758) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(-0.043443505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3869218) q[1];
sx q[1];
rz(-2.6082509) q[1];
sx q[1];
rz(1.5295332) q[1];
x q[2];
rz(2.6838949) q[3];
sx q[3];
rz(-2.5182596) q[3];
sx q[3];
rz(1.8231572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-0.97736248) q[2];
rz(-1.0824925) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.8025773) q[0];
sx q[0];
rz(-0.76887283) q[0];
sx q[0];
rz(-0.73706891) q[0];
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
rz(0.21870273) q[3];
sx q[3];
rz(-2.4730611) q[3];
sx q[3];
rz(-0.685365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
