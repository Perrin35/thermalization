OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35225866) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(2.5369011) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0646677) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(-1.392729) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44126025) q[1];
sx q[1];
rz(-1.2520737) q[1];
sx q[1];
rz(0.031949921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1313685) q[3];
sx q[3];
rz(-1.2636375) q[3];
sx q[3];
rz(2.9205521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.4367746) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(2.2170128) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.3234214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88432089) q[0];
sx q[0];
rz(-1.1182251) q[0];
sx q[0];
rz(1.0898468) q[0];
x q[1];
rz(3.0201868) q[2];
sx q[2];
rz(-0.3519904) q[2];
sx q[2];
rz(2.6999203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9604608) q[1];
sx q[1];
rz(-2.301553) q[1];
sx q[1];
rz(-1.204927) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0074535) q[3];
sx q[3];
rz(-2.4881119) q[3];
sx q[3];
rz(0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(-2.6115131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353697) q[0];
sx q[0];
rz(-0.3100937) q[0];
sx q[0];
rz(1.2139411) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90383756) q[2];
sx q[2];
rz(-2.505216) q[2];
sx q[2];
rz(2.9449376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59280076) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(-0.42216502) q[1];
rz(-pi) q[2];
rz(1.3313053) q[3];
sx q[3];
rz(-1.4756087) q[3];
sx q[3];
rz(2.5483607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(-2.3572671) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.559583) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(2.9942306) q[0];
x q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(-2.6505016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1130484) q[1];
sx q[1];
rz(-1.5469157) q[1];
sx q[1];
rz(-1.6391812) q[1];
rz(1.4284381) q[3];
sx q[3];
rz(-2.5609303) q[3];
sx q[3];
rz(2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-0.56048918) q[2];
rz(3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(-1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1836023) q[0];
sx q[0];
rz(-1.9200268) q[0];
sx q[0];
rz(0.41718418) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1528835) q[2];
sx q[2];
rz(-1.8277797) q[2];
sx q[2];
rz(-0.35494057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69406063) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(1.8156169) q[1];
rz(1.9305265) q[3];
sx q[3];
rz(-1.3621646) q[3];
sx q[3];
rz(-1.9344575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-0.042479854) q[2];
rz(-0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(0.49155864) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(2.65843) q[0];
rz(1.0522316) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(-2.5767456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72950596) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(1.3077523) q[0];
rz(-1.1499314) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.3644621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7397592) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(-3.0763807) q[1];
rz(-pi) q[2];
rz(0.23128831) q[3];
sx q[3];
rz(-1.7219208) q[3];
sx q[3];
rz(-2.1895727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.338039) q[2];
rz(-1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(-0.54779732) q[0];
rz(2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(2.8731667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4329322) q[0];
sx q[0];
rz(-1.4888568) q[0];
sx q[0];
rz(-2.8952778) q[0];
rz(-pi) q[1];
rz(-0.78598778) q[2];
sx q[2];
rz(-2.18835) q[2];
sx q[2];
rz(2.3812889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0722326) q[1];
sx q[1];
rz(-1.4727122) q[1];
sx q[1];
rz(-0.98274883) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8397452) q[3];
sx q[3];
rz(-2.0904623) q[3];
sx q[3];
rz(-2.2083851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(0.8141554) q[2];
rz(0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(2.6928435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1837346) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(0.26259043) q[0];
rz(-pi) q[1];
rz(0.86782311) q[2];
sx q[2];
rz(-0.38735577) q[2];
sx q[2];
rz(-0.58434904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6905578) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(-1.8655538) q[1];
rz(-2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(2.2023647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(0.82434404) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(-0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3138872) q[0];
sx q[0];
rz(-2.8773617) q[0];
sx q[0];
rz(3.0731191) q[0];
rz(-1.47255) q[2];
sx q[2];
rz(-1.7260523) q[2];
sx q[2];
rz(-1.1467903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1259809) q[1];
sx q[1];
rz(-0.19128448) q[1];
sx q[1];
rz(-0.26856883) q[1];
rz(-pi) q[2];
rz(-2.8519565) q[3];
sx q[3];
rz(-1.3445065) q[3];
sx q[3];
rz(2.3895404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5921322) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(-0.67031676) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877514) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(2.642139) q[0];
rz(1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(1.0826899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2912746) q[0];
sx q[0];
rz(-0.58915888) q[0];
sx q[0];
rz(-3.0535112) q[0];
rz(-pi) q[1];
rz(-1.535365) q[2];
sx q[2];
rz(-2.3897768) q[2];
sx q[2];
rz(-2.5978136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.06304) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(-2.5979554) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70268481) q[3];
sx q[3];
rz(-0.48326884) q[3];
sx q[3];
rz(-1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(0.74116771) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(2.1096061) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(0.18556553) q[3];
sx q[3];
rz(-0.75510988) q[3];
sx q[3];
rz(1.818944) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];