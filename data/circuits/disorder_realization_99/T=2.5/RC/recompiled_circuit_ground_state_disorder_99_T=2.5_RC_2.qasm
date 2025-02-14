OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7991601) q[0];
sx q[0];
rz(-2.7437796) q[0];
sx q[0];
rz(1.2678658) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(1.9222577) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3800664) q[0];
sx q[0];
rz(-2.5396569) q[0];
sx q[0];
rz(-1.2720889) q[0];
rz(-pi) q[1];
rz(2.8706495) q[2];
sx q[2];
rz(-0.56519714) q[2];
sx q[2];
rz(-0.2172367) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75324149) q[1];
sx q[1];
rz(-1.2920391) q[1];
sx q[1];
rz(-2.9294347) q[1];
rz(2.9079535) q[3];
sx q[3];
rz(-2.2147182) q[3];
sx q[3];
rz(1.131191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5775602) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(-0.91156256) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0381222) q[0];
sx q[0];
rz(-2.7439674) q[0];
sx q[0];
rz(-3.0225515) q[0];
rz(-1.1301522) q[1];
sx q[1];
rz(-0.96769133) q[1];
sx q[1];
rz(-1.7134604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4864246) q[0];
sx q[0];
rz(-2.0858686) q[0];
sx q[0];
rz(1.0977655) q[0];
x q[1];
rz(2.7409254) q[2];
sx q[2];
rz(-1.9611437) q[2];
sx q[2];
rz(-1.1125178) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94257689) q[1];
sx q[1];
rz(-1.4369094) q[1];
sx q[1];
rz(1.1034637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30562206) q[3];
sx q[3];
rz(-1.8942617) q[3];
sx q[3];
rz(-2.4572069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6156562) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(2.3770135) q[2];
rz(2.6584451) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(-0.84883261) q[3];
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
rz(2.0371542) q[0];
sx q[0];
rz(-2.8762682) q[0];
sx q[0];
rz(1.4720488) q[0];
rz(2.5813591) q[1];
sx q[1];
rz(-1.7543703) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.154523) q[0];
sx q[0];
rz(-1.8049585) q[0];
sx q[0];
rz(-2.9375141) q[0];
rz(-2.4666489) q[2];
sx q[2];
rz(-1.5752572) q[2];
sx q[2];
rz(0.5677815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83207291) q[1];
sx q[1];
rz(-0.64180556) q[1];
sx q[1];
rz(2.1596396) q[1];
rz(-pi) q[2];
rz(1.7974924) q[3];
sx q[3];
rz(-1.8009225) q[3];
sx q[3];
rz(-0.067148681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6560087) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(2.7991925) q[2];
rz(-1.418669) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(-0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(2.2320084) q[0];
sx q[0];
rz(-0.37264687) q[0];
sx q[0];
rz(1.6147856) q[0];
rz(0.80742637) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(-1.2015013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75185173) q[0];
sx q[0];
rz(-0.83503113) q[0];
sx q[0];
rz(2.7324008) q[0];
rz(-pi) q[1];
rz(-0.64862675) q[2];
sx q[2];
rz(-1.0264215) q[2];
sx q[2];
rz(2.8971162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8462584) q[1];
sx q[1];
rz(-2.3148119) q[1];
sx q[1];
rz(2.1200402) q[1];
rz(0.25629429) q[3];
sx q[3];
rz(-1.5206889) q[3];
sx q[3];
rz(-2.9378892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80369192) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(-1.2376415) q[2];
rz(-1.3112274) q[3];
sx q[3];
rz(-1.5179736) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.7685522) q[0];
sx q[0];
rz(-2.232724) q[0];
rz(-0.88242775) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(2.4095101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1279739) q[0];
sx q[0];
rz(-0.92500702) q[0];
sx q[0];
rz(0.60198204) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8516099) q[2];
sx q[2];
rz(-0.78943832) q[2];
sx q[2];
rz(-2.0470999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4914507) q[1];
sx q[1];
rz(-0.60381266) q[1];
sx q[1];
rz(-2.4039388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3957) q[3];
sx q[3];
rz(-0.98222662) q[3];
sx q[3];
rz(2.4323127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.070179209) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(-1.8522235) q[2];
rz(1.6728801) q[3];
sx q[3];
rz(-1.2082992) q[3];
sx q[3];
rz(-2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656723) q[0];
sx q[0];
rz(-1.1812295) q[0];
sx q[0];
rz(1.1015724) q[0];
rz(-2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-2.1563931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6533642) q[0];
sx q[0];
rz(-2.4765827) q[0];
sx q[0];
rz(1.2763766) q[0];
rz(1.9198138) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(1.3427918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8781153) q[1];
sx q[1];
rz(-0.92258555) q[1];
sx q[1];
rz(0.99261673) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5959625) q[3];
sx q[3];
rz(-1.0428793) q[3];
sx q[3];
rz(-2.2759469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78099403) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(-1.194427) q[2];
rz(-3.0418975) q[3];
sx q[3];
rz(-1.7710779) q[3];
sx q[3];
rz(0.97894198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128971) q[0];
sx q[0];
rz(-2.6850057) q[0];
sx q[0];
rz(2.0275443) q[0];
rz(-1.7440589) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(0.31375113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.926696) q[0];
sx q[0];
rz(-0.29698661) q[0];
sx q[0];
rz(1.7420705) q[0];
rz(-pi) q[1];
rz(-2.9767562) q[2];
sx q[2];
rz(-1.4408752) q[2];
sx q[2];
rz(0.31229737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75668979) q[1];
sx q[1];
rz(-2.3339565) q[1];
sx q[1];
rz(-0.4203504) q[1];
rz(-2.8242495) q[3];
sx q[3];
rz(-0.71094027) q[3];
sx q[3];
rz(-1.2430354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.067001192) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(2.5355549) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(-1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17230497) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(1.1267927) q[0];
rz(0.91398319) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(0.21805683) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870616) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(2.0824672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3534691) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(0.1696378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7331024) q[1];
sx q[1];
rz(-2.4839851) q[1];
sx q[1];
rz(1.7773184) q[1];
x q[2];
rz(0.86450926) q[3];
sx q[3];
rz(-2.9614523) q[3];
sx q[3];
rz(0.064398191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(-2.8524354) q[2];
rz(-0.9946) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(-0.4579671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-1.0884322) q[0];
sx q[0];
rz(-0.76643884) q[0];
rz(-1.6038766) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(-0.91805735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1963475) q[0];
sx q[0];
rz(-2.1736988) q[0];
sx q[0];
rz(-0.57749282) q[0];
rz(-1.7154022) q[2];
sx q[2];
rz(-1.40398) q[2];
sx q[2];
rz(-1.5125991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80528211) q[1];
sx q[1];
rz(-2.2754221) q[1];
sx q[1];
rz(-1.8971838) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3206743) q[3];
sx q[3];
rz(-1.6733352) q[3];
sx q[3];
rz(-0.96270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(-0.85203552) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(-2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8155415) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(-0.091212243) q[0];
rz(2.3444029) q[1];
sx q[1];
rz(-1.3858831) q[1];
sx q[1];
rz(0.43112722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6043458) q[0];
sx q[0];
rz(-1.6291999) q[0];
sx q[0];
rz(-1.6184774) q[0];
rz(-pi) q[1];
rz(-0.098478949) q[2];
sx q[2];
rz(-2.2171671) q[2];
sx q[2];
rz(-1.3502094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9111002) q[1];
sx q[1];
rz(-2.6395032) q[1];
sx q[1];
rz(2.6761495) q[1];
x q[2];
rz(-2.1717908) q[3];
sx q[3];
rz(-1.2088115) q[3];
sx q[3];
rz(-2.8014744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5690696) q[2];
sx q[2];
rz(-1.5745682) q[2];
sx q[2];
rz(1.265906) q[2];
rz(-1.8484533) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(0.40614793) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311998) q[0];
sx q[0];
rz(-2.4507903) q[0];
sx q[0];
rz(-2.3660085) q[0];
rz(-0.56397437) q[1];
sx q[1];
rz(-1.0697983) q[1];
sx q[1];
rz(-1.0615798) q[1];
rz(1.3971267) q[2];
sx q[2];
rz(-0.43890719) q[2];
sx q[2];
rz(-0.70499805) q[2];
rz(-2.1544477) q[3];
sx q[3];
rz(-0.55716438) q[3];
sx q[3];
rz(-1.1732994) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
