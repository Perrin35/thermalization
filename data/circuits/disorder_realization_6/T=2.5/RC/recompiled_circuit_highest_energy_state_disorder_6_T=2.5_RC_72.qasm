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
rz(-1.8974798) q[0];
sx q[0];
rz(-0.31960684) q[0];
sx q[0];
rz(2.6407916) q[0];
rz(-2.8473941) q[1];
sx q[1];
rz(-0.74217141) q[1];
sx q[1];
rz(-0.54744005) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1337903) q[0];
sx q[0];
rz(-0.34388516) q[0];
sx q[0];
rz(2.6965428) q[0];
x q[1];
rz(2.6771678) q[2];
sx q[2];
rz(-1.8193442) q[2];
sx q[2];
rz(-0.82962245) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3708122) q[1];
sx q[1];
rz(-2.6264084) q[1];
sx q[1];
rz(2.9574802) q[1];
rz(-pi) q[2];
rz(2.1294598) q[3];
sx q[3];
rz(-0.61311537) q[3];
sx q[3];
rz(2.0215542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1529634) q[2];
sx q[2];
rz(-0.84501481) q[2];
sx q[2];
rz(0.60626924) q[2];
rz(-1.2074977) q[3];
sx q[3];
rz(-1.8469801) q[3];
sx q[3];
rz(1.5509563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2400874) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(-2.4497581) q[0];
rz(1.8531307) q[1];
sx q[1];
rz(-1.144578) q[1];
sx q[1];
rz(0.28786927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3436056) q[0];
sx q[0];
rz(-1.4427066) q[0];
sx q[0];
rz(-1.8592181) q[0];
rz(-2.1580196) q[2];
sx q[2];
rz(-2.1701239) q[2];
sx q[2];
rz(0.22967185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.560656) q[1];
sx q[1];
rz(-2.0163149) q[1];
sx q[1];
rz(-1.5195283) q[1];
rz(-0.13707844) q[3];
sx q[3];
rz(-1.6774639) q[3];
sx q[3];
rz(-0.24124537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3466779) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(2.1686926) q[2];
rz(1.8074624) q[3];
sx q[3];
rz(-0.79136807) q[3];
sx q[3];
rz(-2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7767104) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(-2.8926988) q[0];
rz(-1.6913255) q[1];
sx q[1];
rz(-1.1509117) q[1];
sx q[1];
rz(0.90363622) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2657651) q[0];
sx q[0];
rz(-0.96185124) q[0];
sx q[0];
rz(-1.5399571) q[0];
rz(-pi) q[1];
rz(-1.9157639) q[2];
sx q[2];
rz(-1.4699226) q[2];
sx q[2];
rz(-0.28988923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94374471) q[1];
sx q[1];
rz(-1.1863803) q[1];
sx q[1];
rz(0.77898394) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60402212) q[3];
sx q[3];
rz(-1.2327177) q[3];
sx q[3];
rz(0.018751831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9768208) q[2];
sx q[2];
rz(-0.73126078) q[2];
sx q[2];
rz(2.8110647) q[2];
rz(-0.22322379) q[3];
sx q[3];
rz(-2.0285172) q[3];
sx q[3];
rz(-2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6966003) q[0];
sx q[0];
rz(-1.3560504) q[0];
sx q[0];
rz(0.14183216) q[0];
rz(0.089135535) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(-0.95318046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2640799) q[0];
sx q[0];
rz(-0.82897128) q[0];
sx q[0];
rz(2.6116051) q[0];
x q[1];
rz(2.8492979) q[2];
sx q[2];
rz(-0.39543786) q[2];
sx q[2];
rz(1.866445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5281727) q[1];
sx q[1];
rz(-2.0866359) q[1];
sx q[1];
rz(-0.9524166) q[1];
x q[2];
rz(2.6450883) q[3];
sx q[3];
rz(-1.3937573) q[3];
sx q[3];
rz(-0.6547218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7715093) q[2];
sx q[2];
rz(-3.0107095) q[2];
sx q[2];
rz(-0.5292325) q[2];
rz(-1.0445163) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(-0.15531003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.8998242) q[0];
sx q[0];
rz(-0.45503765) q[0];
rz(0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(0.98247772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331086) q[0];
sx q[0];
rz(-0.31705515) q[0];
sx q[0];
rz(0.32755537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58042629) q[2];
sx q[2];
rz(-2.1820445) q[2];
sx q[2];
rz(-2.9410604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8807672) q[1];
sx q[1];
rz(-1.2952571) q[1];
sx q[1];
rz(1.6795621) q[1];
x q[2];
rz(-0.69070821) q[3];
sx q[3];
rz(-0.36223199) q[3];
sx q[3];
rz(-1.9286523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5386397) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(0.61720103) q[2];
rz(1.3868825) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(-2.9673747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0844326) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(-0.97849751) q[0];
rz(-2.1242566) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(2.691958) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89116865) q[0];
sx q[0];
rz(-0.40968597) q[0];
sx q[0];
rz(0.72575672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.041018) q[2];
sx q[2];
rz(-0.82056773) q[2];
sx q[2];
rz(-2.5448397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3842135) q[1];
sx q[1];
rz(-1.6156726) q[1];
sx q[1];
rz(1.0991389) q[1];
rz(-pi) q[2];
rz(2.4586304) q[3];
sx q[3];
rz(-2.0724845) q[3];
sx q[3];
rz(1.7428118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6189239) q[2];
sx q[2];
rz(-1.3441244) q[2];
sx q[2];
rz(-1.3783003) q[2];
rz(-1.1528692) q[3];
sx q[3];
rz(-1.9882354) q[3];
sx q[3];
rz(-2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744837) q[0];
sx q[0];
rz(-0.672988) q[0];
sx q[0];
rz(0.33669499) q[0];
rz(0.7252655) q[1];
sx q[1];
rz(-1.5870353) q[1];
sx q[1];
rz(2.7560962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9850901) q[0];
sx q[0];
rz(-0.96961951) q[0];
sx q[0];
rz(1.1416092) q[0];
x q[1];
rz(-2.9988937) q[2];
sx q[2];
rz(-1.8178504) q[2];
sx q[2];
rz(-2.3076535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4591864) q[1];
sx q[1];
rz(-1.2810315) q[1];
sx q[1];
rz(-2.2166445) q[1];
rz(0.5637873) q[3];
sx q[3];
rz(-2.4372589) q[3];
sx q[3];
rz(-0.35455656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5943299) q[2];
sx q[2];
rz(-1.2196093) q[2];
sx q[2];
rz(-0.29898137) q[2];
rz(-0.79214054) q[3];
sx q[3];
rz(-2.7463089) q[3];
sx q[3];
rz(-1.3778752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250273) q[0];
sx q[0];
rz(-1.4341609) q[0];
sx q[0];
rz(-0.7780956) q[0];
rz(-1.7447225) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(-0.088833749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0087826) q[0];
sx q[0];
rz(-2.7799921) q[0];
sx q[0];
rz(2.8426991) q[0];
rz(-2.2065032) q[2];
sx q[2];
rz(-0.82391058) q[2];
sx q[2];
rz(1.3694135) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1792887) q[1];
sx q[1];
rz(-1.3109164) q[1];
sx q[1];
rz(2.0043027) q[1];
rz(-pi) q[2];
rz(-3.1109654) q[3];
sx q[3];
rz(-1.9383241) q[3];
sx q[3];
rz(2.5961329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54625154) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(-2.0383932) q[2];
rz(2.5698419) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(1.5248689) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095846) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(0.65015656) q[0];
rz(2.0434168) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(0.065902725) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377334) q[0];
sx q[0];
rz(-2.2324076) q[0];
sx q[0];
rz(2.1801722) q[0];
rz(2.8464855) q[2];
sx q[2];
rz(-2.9431651) q[2];
sx q[2];
rz(1.0498384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2615589) q[1];
sx q[1];
rz(-1.7249692) q[1];
sx q[1];
rz(-2.3460813) q[1];
x q[2];
rz(1.8239924) q[3];
sx q[3];
rz(-2.2166219) q[3];
sx q[3];
rz(-1.822804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1419475) q[2];
sx q[2];
rz(-1.9530752) q[2];
sx q[2];
rz(-0.56335062) q[2];
rz(0.79936409) q[3];
sx q[3];
rz(-2.1269709) q[3];
sx q[3];
rz(2.8857901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946741) q[0];
sx q[0];
rz(-3.1212555) q[0];
sx q[0];
rz(2.643423) q[0];
rz(0.26214552) q[1];
sx q[1];
rz(-1.7356977) q[1];
sx q[1];
rz(1.3462876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57980832) q[0];
sx q[0];
rz(-2.0286848) q[0];
sx q[0];
rz(-2.3696128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.002368517) q[2];
sx q[2];
rz(-1.5813229) q[2];
sx q[2];
rz(-2.0879371) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7340103) q[1];
sx q[1];
rz(-2.2395113) q[1];
sx q[1];
rz(2.0120502) q[1];
rz(-pi) q[2];
rz(2.3614083) q[3];
sx q[3];
rz(-1.4139851) q[3];
sx q[3];
rz(-3.1136857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6654309) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(-2.2594182) q[2];
rz(-0.65794182) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(2.1420124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6366202) q[0];
sx q[0];
rz(-1.636878) q[0];
sx q[0];
rz(-1.0374877) q[0];
rz(2.6344015) q[1];
sx q[1];
rz(-2.3585944) q[1];
sx q[1];
rz(2.0814887) q[1];
rz(0.29303356) q[2];
sx q[2];
rz(-1.3923981) q[2];
sx q[2];
rz(-0.57033718) q[2];
rz(-0.67340452) q[3];
sx q[3];
rz(-0.87457228) q[3];
sx q[3];
rz(1.2403929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
