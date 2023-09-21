OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.7339098) q[0];
sx q[0];
rz(1.9440584) q[0];
rz(-pi) q[1];
rz(1.5288562) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(0.97181335) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.672294) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(-1.6507571) q[1];
rz(1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(-1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97074189) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(3.0424776) q[0];
rz(-0.72950659) q[2];
sx q[2];
rz(-1.7765877) q[2];
sx q[2];
rz(-1.4480928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.193589) q[1];
sx q[1];
rz(-1.6063599) q[1];
sx q[1];
rz(2.8405632) q[1];
rz(-0.64579441) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16386579) q[0];
sx q[0];
rz(-2.8054872) q[0];
sx q[0];
rz(-1.1019812) q[0];
x q[1];
rz(0.83061647) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-0.22731552) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1516583) q[1];
sx q[1];
rz(-0.5608359) q[1];
sx q[1];
rz(2.6480617) q[1];
x q[2];
rz(-2.0796892) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(-3.0825465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(-2.3245658) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(0.27483637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66514689) q[0];
sx q[0];
rz(-0.11867141) q[0];
sx q[0];
rz(2.2975886) q[0];
rz(3.1399973) q[2];
sx q[2];
rz(-1.110382) q[2];
sx q[2];
rz(2.7574725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8400164) q[1];
sx q[1];
rz(-0.77743545) q[1];
sx q[1];
rz(2.2393054) q[1];
rz(0.87423012) q[3];
sx q[3];
rz(-1.6682373) q[3];
sx q[3];
rz(-1.3144573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65790025) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.6246187) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.330634) q[0];
sx q[0];
rz(-1.2067544) q[0];
sx q[0];
rz(2.2383658) q[0];
rz(-0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.595572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-0.43069559) q[1];
rz(0.77002854) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-2.8175763) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.5312622) q[3];
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
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(-2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-0.34067672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9522889) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(-3.1097263) q[0];
rz(-pi) q[1];
rz(-1.864205) q[2];
sx q[2];
rz(-2.3050606) q[2];
sx q[2];
rz(-2.4794527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50960474) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(-1.3118841) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42231456) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(-2.5781393) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(-0.46494928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5349605) q[0];
sx q[0];
rz(-0.50919845) q[0];
sx q[0];
rz(1.9726994) q[0];
rz(-pi) q[1];
rz(0.99636997) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(-0.043957274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9629434) q[1];
sx q[1];
rz(-1.4649676) q[1];
sx q[1];
rz(-2.9220198) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0602337) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(0.57146227) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(3.0632339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6358444) q[0];
sx q[0];
rz(-0.76857476) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80348357) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(1.6840881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9662712) q[1];
sx q[1];
rz(-2.2178855) q[1];
sx q[1];
rz(-0.21827571) q[1];
x q[2];
rz(1.4140698) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(-2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(0.87402469) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6459991) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(-1.5828703) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36058493) q[2];
sx q[2];
rz(-1.9546095) q[2];
sx q[2];
rz(-1.4521445) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2487138) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(-1.0700657) q[1];
rz(-0.47253982) q[3];
sx q[3];
rz(-2.4774385) q[3];
sx q[3];
rz(3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.727227) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(-1.1788517) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154685) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(-1.8490851) q[0];
rz(-pi) q[1];
rz(1.6640501) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(1.7699514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8995754) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(2.1980397) q[1];
rz(-pi) q[2];
rz(-2.6188649) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(-2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.4609059) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.4964676) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(1.7189797) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
