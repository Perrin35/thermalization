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
rz(1.2441128) q[0];
sx q[0];
rz(-2.8219858) q[0];
sx q[0];
rz(0.50080103) q[0];
rz(0.29419857) q[1];
sx q[1];
rz(3.8837641) q[1];
sx q[1];
rz(9.972218) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1337903) q[0];
sx q[0];
rz(-2.7977075) q[0];
sx q[0];
rz(-2.6965428) q[0];
rz(-pi) q[1];
rz(0.4644248) q[2];
sx q[2];
rz(-1.8193442) q[2];
sx q[2];
rz(-2.3119702) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55996229) q[1];
sx q[1];
rz(-2.0764372) q[1];
sx q[1];
rz(1.4675114) q[1];
rz(-pi) q[2];
rz(1.0121328) q[3];
sx q[3];
rz(-0.61311537) q[3];
sx q[3];
rz(1.1200384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1529634) q[2];
sx q[2];
rz(-2.2965778) q[2];
sx q[2];
rz(-0.60626924) q[2];
rz(-1.9340949) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(-1.5906364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2400874) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(-0.69183451) q[0];
rz(-1.2884619) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(-0.28786927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.797987) q[0];
sx q[0];
rz(-1.6988861) q[0];
sx q[0];
rz(1.8592181) q[0];
rz(-pi) q[1];
rz(0.68142724) q[2];
sx q[2];
rz(-2.3286901) q[2];
sx q[2];
rz(2.5035448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1296245) q[1];
sx q[1];
rz(-1.6170562) q[1];
sx q[1];
rz(-2.6955626) q[1];
x q[2];
rz(-3.0045142) q[3];
sx q[3];
rz(-1.4641288) q[3];
sx q[3];
rz(-0.24124537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3466779) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(-2.1686926) q[2];
rz(-1.3341303) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(-0.62132519) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36488229) q[0];
sx q[0];
rz(-0.24672395) q[0];
sx q[0];
rz(-0.2488939) q[0];
rz(-1.4502672) q[1];
sx q[1];
rz(-1.1509117) q[1];
sx q[1];
rz(2.2379564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2657651) q[0];
sx q[0];
rz(-2.1797414) q[0];
sx q[0];
rz(1.5399571) q[0];
rz(-pi) q[1];
rz(0.10714178) q[2];
sx q[2];
rz(-1.2276548) q[2];
sx q[2];
rz(1.3170774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26402707) q[1];
sx q[1];
rz(-2.2912044) q[1];
sx q[1];
rz(-2.6191465) q[1];
x q[2];
rz(-2.5375705) q[3];
sx q[3];
rz(-1.2327177) q[3];
sx q[3];
rz(3.1228408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9768208) q[2];
sx q[2];
rz(-2.4103319) q[2];
sx q[2];
rz(-0.33052793) q[2];
rz(0.22322379) q[3];
sx q[3];
rz(-1.1130755) q[3];
sx q[3];
rz(-2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44499236) q[0];
sx q[0];
rz(-1.3560504) q[0];
sx q[0];
rz(2.9997605) q[0];
rz(-0.089135535) q[1];
sx q[1];
rz(-2.8991957) q[1];
sx q[1];
rz(-0.95318046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9784666) q[0];
sx q[0];
rz(-0.88136601) q[0];
sx q[0];
rz(1.0667144) q[0];
x q[1];
rz(1.6905027) q[2];
sx q[2];
rz(-1.1929995) q[2];
sx q[2];
rz(2.1816513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.492748) q[1];
sx q[1];
rz(-2.358583) q[1];
sx q[1];
rz(-0.79644189) q[1];
rz(-pi) q[2];
rz(2.7823086) q[3];
sx q[3];
rz(-0.52462884) q[3];
sx q[3];
rz(0.60177647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7715093) q[2];
sx q[2];
rz(-3.0107095) q[2];
sx q[2];
rz(0.5292325) q[2];
rz(-1.0445163) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(-0.15531003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.2417685) q[0];
sx q[0];
rz(-2.686555) q[0];
rz(0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(0.98247772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7916103) q[0];
sx q[0];
rz(-1.6712708) q[0];
sx q[0];
rz(-2.8403681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90679866) q[2];
sx q[2];
rz(-2.325146) q[2];
sx q[2];
rz(-2.0891669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12085401) q[1];
sx q[1];
rz(-2.845872) q[1];
sx q[1];
rz(-2.775008) q[1];
rz(1.3339046) q[3];
sx q[3];
rz(-1.8474527) q[3];
sx q[3];
rz(1.93678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.602953) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(2.5243916) q[2];
rz(1.3868825) q[3];
sx q[3];
rz(-2.230547) q[3];
sx q[3];
rz(-0.174218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0844326) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(-2.1630951) q[0];
rz(2.1242566) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(0.4496347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7788197) q[0];
sx q[0];
rz(-1.8383433) q[0];
sx q[0];
rz(-0.31407251) q[0];
x q[1];
rz(1.1005747) q[2];
sx q[2];
rz(-2.3210249) q[2];
sx q[2];
rz(-0.596753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3510531) q[1];
sx q[1];
rz(-2.0419403) q[1];
sx q[1];
rz(-3.0912249) q[1];
rz(-pi) q[2];
rz(-2.4261023) q[3];
sx q[3];
rz(-0.82287298) q[3];
sx q[3];
rz(0.70586849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6189239) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(-1.3783003) q[2];
rz(1.9887234) q[3];
sx q[3];
rz(-1.1533573) q[3];
sx q[3];
rz(-0.26067477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744837) q[0];
sx q[0];
rz(-0.672988) q[0];
sx q[0];
rz(0.33669499) q[0];
rz(-2.4163272) q[1];
sx q[1];
rz(-1.5870353) q[1];
sx q[1];
rz(-0.38549647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3025997) q[0];
sx q[0];
rz(-1.9210879) q[0];
sx q[0];
rz(2.4953609) q[0];
rz(1.0573559) q[2];
sx q[2];
rz(-0.2845736) q[2];
sx q[2];
rz(1.7764661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8908317) q[1];
sx q[1];
rz(-2.4422993) q[1];
sx q[1];
rz(-1.1108524) q[1];
rz(-pi) q[2];
rz(1.1445505) q[3];
sx q[3];
rz(-0.9916456) q[3];
sx q[3];
rz(-0.338011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.54726279) q[2];
sx q[2];
rz(-1.2196093) q[2];
sx q[2];
rz(2.8426113) q[2];
rz(0.79214054) q[3];
sx q[3];
rz(-0.39528379) q[3];
sx q[3];
rz(-1.3778752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91656536) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(2.3634971) q[0];
rz(1.7447225) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(-3.0527589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13281001) q[0];
sx q[0];
rz(-2.7799921) q[0];
sx q[0];
rz(0.2988935) q[0];
x q[1];
rz(-0.93508945) q[2];
sx q[2];
rz(-0.82391058) q[2];
sx q[2];
rz(-1.3694135) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11528291) q[1];
sx q[1];
rz(-0.50116083) q[1];
sx q[1];
rz(1.0064685) q[1];
rz(-pi) q[2];
rz(1.2031112) q[3];
sx q[3];
rz(-1.5993777) q[3];
sx q[3];
rz(-2.1052484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5953411) q[2];
sx q[2];
rz(-2.310414) q[2];
sx q[2];
rz(1.1031995) q[2];
rz(-0.57175076) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(-1.6167238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53200805) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(-2.4914361) q[0];
rz(2.0434168) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(0.065902725) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1797721) q[0];
sx q[0];
rz(-2.0393436) q[0];
sx q[0];
rz(-0.75956068) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9515334) q[2];
sx q[2];
rz(-1.628161) q[2];
sx q[2];
rz(0.8105958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3014631) q[1];
sx q[1];
rz(-2.3345229) q[1];
sx q[1];
rz(0.21424861) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4801257) q[3];
sx q[3];
rz(-1.3693878) q[3];
sx q[3];
rz(-0.097518343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1419475) q[2];
sx q[2];
rz(-1.1885175) q[2];
sx q[2];
rz(-0.56335062) q[2];
rz(0.79936409) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(0.25580251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946741) q[0];
sx q[0];
rz(-3.1212555) q[0];
sx q[0];
rz(-0.49816966) q[0];
rz(2.8794471) q[1];
sx q[1];
rz(-1.7356977) q[1];
sx q[1];
rz(1.795305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243626) q[0];
sx q[0];
rz(-0.87275332) q[0];
sx q[0];
rz(-2.5265187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5602697) q[2];
sx q[2];
rz(-1.5684279) q[2];
sx q[2];
rz(-2.6244768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.263239) q[1];
sx q[1];
rz(-1.9124971) q[1];
sx q[1];
rz(-0.71818761) q[1];
x q[2];
rz(-0.22110053) q[3];
sx q[3];
rz(-2.3490864) q[3];
sx q[3];
rz(-1.4421876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6654309) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(2.2594182) q[2];
rz(0.65794182) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(0.9995802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5049725) q[0];
sx q[0];
rz(-1.636878) q[0];
sx q[0];
rz(-1.0374877) q[0];
rz(-0.50719117) q[1];
sx q[1];
rz(-2.3585944) q[1];
sx q[1];
rz(2.0814887) q[1];
rz(0.55804923) q[2];
sx q[2];
rz(-0.34172575) q[2];
sx q[2];
rz(1.5321129) q[2];
rz(2.389659) q[3];
sx q[3];
rz(-1.0718413) q[3];
sx q[3];
rz(-2.9990673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
