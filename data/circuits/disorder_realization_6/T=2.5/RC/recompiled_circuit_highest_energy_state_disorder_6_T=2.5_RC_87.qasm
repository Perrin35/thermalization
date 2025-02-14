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
rz(0.29419857) q[1];
sx q[1];
rz(-2.3994212) q[1];
sx q[1];
rz(0.54744005) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1337903) q[0];
sx q[0];
rz(-2.7977075) q[0];
sx q[0];
rz(-0.44504984) q[0];
rz(-pi) q[1];
rz(-0.51551657) q[2];
sx q[2];
rz(-0.52243865) q[2];
sx q[2];
rz(-1.197627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5816304) q[1];
sx q[1];
rz(-1.0651555) q[1];
sx q[1];
rz(1.6740812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7846401) q[3];
sx q[3];
rz(-1.0610749) q[3];
sx q[3];
rz(-1.3689643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1529634) q[2];
sx q[2];
rz(-0.84501481) q[2];
sx q[2];
rz(-0.60626924) q[2];
rz(-1.9340949) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(1.5509563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9015053) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(0.69183451) q[0];
rz(-1.2884619) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(2.8537234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9622274) q[0];
sx q[0];
rz(-0.31485885) q[0];
sx q[0];
rz(-1.1456144) q[0];
x q[1];
rz(-2.4601654) q[2];
sx q[2];
rz(-0.81290258) q[2];
sx q[2];
rz(0.63804783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1296245) q[1];
sx q[1];
rz(-1.5245364) q[1];
sx q[1];
rz(2.6955626) q[1];
rz(-pi) q[2];
rz(1.6784662) q[3];
sx q[3];
rz(-1.4345019) q[3];
sx q[3];
rz(1.8267269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3466779) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(2.1686926) q[2];
rz(-1.8074624) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(0.62132519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36488229) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(-0.2488939) q[0];
rz(1.6913255) q[1];
sx q[1];
rz(-1.9906809) q[1];
sx q[1];
rz(0.90363622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2657651) q[0];
sx q[0];
rz(-2.1797414) q[0];
sx q[0];
rz(1.6016356) q[0];
rz(1.2258287) q[2];
sx q[2];
rz(-1.6716701) q[2];
sx q[2];
rz(0.28988923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.26402707) q[1];
sx q[1];
rz(-0.85038821) q[1];
sx q[1];
rz(2.6191465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9744942) q[3];
sx q[3];
rz(-1.0053325) q[3];
sx q[3];
rz(-1.3270449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9768208) q[2];
sx q[2];
rz(-2.4103319) q[2];
sx q[2];
rz(2.8110647) q[2];
rz(2.9183689) q[3];
sx q[3];
rz(-2.0285172) q[3];
sx q[3];
rz(-2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6966003) q[0];
sx q[0];
rz(-1.7855423) q[0];
sx q[0];
rz(0.14183216) q[0];
rz(3.0524571) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(0.95318046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16312604) q[0];
sx q[0];
rz(-0.88136601) q[0];
sx q[0];
rz(-2.0748782) q[0];
x q[1];
rz(2.8492979) q[2];
sx q[2];
rz(-0.39543786) q[2];
sx q[2];
rz(1.866445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5281727) q[1];
sx q[1];
rz(-1.0549567) q[1];
sx q[1];
rz(-0.9524166) q[1];
x q[2];
rz(0.49650438) q[3];
sx q[3];
rz(-1.7478353) q[3];
sx q[3];
rz(2.4868709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7715093) q[2];
sx q[2];
rz(-3.0107095) q[2];
sx q[2];
rz(2.6123602) q[2];
rz(-1.0445163) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(2.9862826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.2417685) q[0];
sx q[0];
rz(-2.686555) q[0];
rz(-0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(2.1591149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.608484) q[0];
sx q[0];
rz(-0.31705515) q[0];
sx q[0];
rz(-2.8140373) q[0];
rz(-pi) q[1];
rz(2.2682954) q[2];
sx q[2];
rz(-1.1050537) q[2];
sx q[2];
rz(-1.0103152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0207386) q[1];
sx q[1];
rz(-0.29572067) q[1];
sx q[1];
rz(2.775008) q[1];
rz(-1.3339046) q[3];
sx q[3];
rz(-1.8474527) q[3];
sx q[3];
rz(1.2048126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5386397) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(0.61720103) q[2];
rz(-1.3868825) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(-0.174218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0844326) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(2.1630951) q[0];
rz(2.1242566) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(-2.691958) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627729) q[0];
sx q[0];
rz(-1.3032493) q[0];
sx q[0];
rz(-0.31407251) q[0];
rz(-pi) q[1];
rz(-1.1005747) q[2];
sx q[2];
rz(-0.82056773) q[2];
sx q[2];
rz(-0.596753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3510531) q[1];
sx q[1];
rz(-2.0419403) q[1];
sx q[1];
rz(-0.050367722) q[1];
x q[2];
rz(-0.95532598) q[3];
sx q[3];
rz(-2.1571473) q[3];
sx q[3];
rz(-0.20092703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52266878) q[2];
sx q[2];
rz(-1.3441244) q[2];
sx q[2];
rz(-1.3783003) q[2];
rz(1.9887234) q[3];
sx q[3];
rz(-1.9882354) q[3];
sx q[3];
rz(-2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5744837) q[0];
sx q[0];
rz(-2.4686047) q[0];
sx q[0];
rz(0.33669499) q[0];
rz(-2.4163272) q[1];
sx q[1];
rz(-1.5870353) q[1];
sx q[1];
rz(-0.38549647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1565026) q[0];
sx q[0];
rz(-0.96961951) q[0];
sx q[0];
rz(1.1416092) q[0];
rz(-pi) q[1];
rz(-2.9988937) q[2];
sx q[2];
rz(-1.3237423) q[2];
sx q[2];
rz(-0.83393914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8908317) q[1];
sx q[1];
rz(-2.4422993) q[1];
sx q[1];
rz(2.0307402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5637873) q[3];
sx q[3];
rz(-2.4372589) q[3];
sx q[3];
rz(2.7870361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5943299) q[2];
sx q[2];
rz(-1.9219834) q[2];
sx q[2];
rz(0.29898137) q[2];
rz(2.3494521) q[3];
sx q[3];
rz(-0.39528379) q[3];
sx q[3];
rz(-1.7637174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91656536) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(-0.7780956) q[0];
rz(1.3968702) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(3.0527589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0087826) q[0];
sx q[0];
rz(-0.36160053) q[0];
sx q[0];
rz(-2.8426991) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85529144) q[2];
sx q[2];
rz(-1.119985) q[2];
sx q[2];
rz(-0.26329783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6316996) q[1];
sx q[1];
rz(-1.1527779) q[1];
sx q[1];
rz(2.8565744) q[1];
rz(-pi) q[2];
rz(3.1109654) q[3];
sx q[3];
rz(-1.2032685) q[3];
sx q[3];
rz(2.5961329) q[3];
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
rz(-2.5751556) q[3];
sx q[3];
rz(1.6167238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095846) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(0.65015656) q[0];
rz(-2.0434168) q[1];
sx q[1];
rz(-1.2547341) q[1];
sx q[1];
rz(0.065902725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1797721) q[0];
sx q[0];
rz(-1.1022491) q[0];
sx q[0];
rz(-2.382032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6292105) q[2];
sx q[2];
rz(-1.3810535) q[2];
sx q[2];
rz(-2.3924215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3014631) q[1];
sx q[1];
rz(-0.80706978) q[1];
sx q[1];
rz(0.21424861) q[1];
rz(-pi) q[2];
rz(-1.8239924) q[3];
sx q[3];
rz(-0.92497073) q[3];
sx q[3];
rz(1.3187887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1419475) q[2];
sx q[2];
rz(-1.1885175) q[2];
sx q[2];
rz(2.578242) q[2];
rz(-2.3422286) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(0.25580251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19485165) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(2.643423) q[0];
rz(-0.26214552) q[1];
sx q[1];
rz(-1.7356977) q[1];
sx q[1];
rz(-1.3462876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57980832) q[0];
sx q[0];
rz(-2.0286848) q[0];
sx q[0];
rz(2.3696128) q[0];
rz(-0.002368517) q[2];
sx q[2];
rz(-1.5602698) q[2];
sx q[2];
rz(-2.0879371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0829858) q[1];
sx q[1];
rz(-2.3595361) q[1];
sx q[1];
rz(2.6460827) q[1];
rz(-pi) q[2];
rz(-0.22110053) q[3];
sx q[3];
rz(-0.79250626) q[3];
sx q[3];
rz(1.4421876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4761618) q[2];
sx q[2];
rz(-2.4527145) q[2];
sx q[2];
rz(-0.8821744) q[2];
rz(-0.65794182) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(2.1420124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(2.6344015) q[1];
sx q[1];
rz(-2.3585944) q[1];
sx q[1];
rz(2.0814887) q[1];
rz(-0.29303356) q[2];
sx q[2];
rz(-1.7491946) q[2];
sx q[2];
rz(2.5712555) q[2];
rz(0.75193361) q[3];
sx q[3];
rz(-2.0697513) q[3];
sx q[3];
rz(0.14252539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
