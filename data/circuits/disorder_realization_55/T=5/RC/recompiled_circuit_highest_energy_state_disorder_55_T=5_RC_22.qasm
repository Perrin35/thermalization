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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(-2.1314148) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(-2.8316166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77783647) q[0];
sx q[0];
rz(-1.68206) q[0];
sx q[0];
rz(2.7194752) q[0];
x q[1];
rz(-1.5484137) q[2];
sx q[2];
rz(-1.7760008) q[2];
sx q[2];
rz(2.4930645) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3821564) q[1];
sx q[1];
rz(-1.8983937) q[1];
sx q[1];
rz(-2.9345678) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83155379) q[3];
sx q[3];
rz(-1.420701) q[3];
sx q[3];
rz(0.41285601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6628722) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-3.0459246) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(-1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9005168) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(-0.61479968) q[0];
rz(0.88848937) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(2.6420171) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1785806) q[0];
sx q[0];
rz(-2.7546429) q[0];
sx q[0];
rz(2.1301756) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4128628) q[2];
sx q[2];
rz(-1.0519895) q[2];
sx q[2];
rz(-0.11473303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.080729) q[1];
sx q[1];
rz(-1.398096) q[1];
sx q[1];
rz(0.91367803) q[1];
x q[2];
rz(2.0604443) q[3];
sx q[3];
rz(-0.64125618) q[3];
sx q[3];
rz(-1.7288417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2483612) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(0.020817967) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(1.9564691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(-0.33367208) q[0];
rz(0.73792136) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(-0.12942448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564011) q[0];
sx q[0];
rz(-1.7298152) q[0];
sx q[0];
rz(-0.88944737) q[0];
rz(-2.470181) q[2];
sx q[2];
rz(-1.08403) q[2];
sx q[2];
rz(-2.0109796) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94679615) q[1];
sx q[1];
rz(-1.692564) q[1];
sx q[1];
rz(-3.0508079) q[1];
rz(1.09472) q[3];
sx q[3];
rz(-2.2942846) q[3];
sx q[3];
rz(-2.0499961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(1.7714436) q[2];
rz(-2.395199) q[3];
sx q[3];
rz(-1.3179702) q[3];
sx q[3];
rz(-0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(-0.56513894) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(-0.82829222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1734742) q[0];
sx q[0];
rz(-0.83195639) q[0];
sx q[0];
rz(-1.678874) q[0];
rz(-2.8548334) q[2];
sx q[2];
rz(-2.0398524) q[2];
sx q[2];
rz(-0.048407528) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2678623) q[1];
sx q[1];
rz(-2.5791188) q[1];
sx q[1];
rz(2.6270694) q[1];
rz(-pi) q[2];
rz(-0.24803646) q[3];
sx q[3];
rz(-1.0674849) q[3];
sx q[3];
rz(-2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68507489) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(-1.3294719) q[0];
rz(-1.9895408) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.7547742) q[0];
sx q[0];
rz(-1.3217682) q[0];
rz(-1.0365965) q[2];
sx q[2];
rz(-1.4956258) q[2];
sx q[2];
rz(0.49803621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21890404) q[1];
sx q[1];
rz(-0.62059015) q[1];
sx q[1];
rz(-0.94938486) q[1];
rz(3.0609691) q[3];
sx q[3];
rz(-2.4084457) q[3];
sx q[3];
rz(-2.6997363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0033215) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(-3.1357583) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(1.4877315) q[0];
rz(-0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(-0.94246513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2104032) q[0];
sx q[0];
rz(-1.382834) q[0];
sx q[0];
rz(0.24954777) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8428188) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(0.58591671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8854495) q[1];
sx q[1];
rz(-2.0383308) q[1];
sx q[1];
rz(-0.3216775) q[1];
x q[2];
rz(0.26047892) q[3];
sx q[3];
rz(-2.0325869) q[3];
sx q[3];
rz(2.8094629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59948644) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(-1.8360809) q[2];
rz(-2.5944338) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9611573) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(0.10051522) q[0];
rz(-0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(0.85711342) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2861569) q[0];
sx q[0];
rz(-1.248994) q[0];
sx q[0];
rz(0.67097539) q[0];
rz(-0.64684644) q[2];
sx q[2];
rz(-1.610272) q[2];
sx q[2];
rz(-1.0919065) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3639314) q[1];
sx q[1];
rz(-2.6865733) q[1];
sx q[1];
rz(-2.6520686) q[1];
rz(0.25730142) q[3];
sx q[3];
rz(-1.0142148) q[3];
sx q[3];
rz(2.0428994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(2.1619469) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7804467) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(3.1106023) q[0];
rz(-0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.8642289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5073858) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(-0.35315634) q[0];
x q[1];
rz(-1.1436815) q[2];
sx q[2];
rz(-2.6027205) q[2];
sx q[2];
rz(-0.32921916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7614224) q[1];
sx q[1];
rz(-1.644165) q[1];
sx q[1];
rz(-2.0580893) q[1];
x q[2];
rz(0.061789024) q[3];
sx q[3];
rz(-1.8553858) q[3];
sx q[3];
rz(0.20670836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.10451) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697407) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(2.790614) q[0];
rz(-0.17414302) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(1.4016271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021371776) q[0];
sx q[0];
rz(-1.8184156) q[0];
sx q[0];
rz(1.4445452) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0171595) q[2];
sx q[2];
rz(-0.87136641) q[2];
sx q[2];
rz(1.0126737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8216711) q[1];
sx q[1];
rz(-0.18583365) q[1];
sx q[1];
rz(1.4129078) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74123592) q[3];
sx q[3];
rz(-2.0307856) q[3];
sx q[3];
rz(-2.7335087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-0.6366716) q[2];
rz(-1.1104256) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(-2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670032) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(-1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.6148753) q[1];
sx q[1];
rz(-2.0775332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23180873) q[0];
sx q[0];
rz(-1.5752821) q[0];
sx q[0];
rz(1.5741328) q[0];
rz(-2.3346155) q[2];
sx q[2];
rz(-0.39013559) q[2];
sx q[2];
rz(1.5240508) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1326931) q[1];
sx q[1];
rz(-1.510396) q[1];
sx q[1];
rz(1.5302883) q[1];
rz(3.0711932) q[3];
sx q[3];
rz(-0.46050763) q[3];
sx q[3];
rz(-1.2706336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4260063) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(3.0675724) q[2];
rz(-2.7600539) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1114125) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(-1.9520957) q[2];
sx q[2];
rz(-1.0972037) q[2];
sx q[2];
rz(0.9158132) q[2];
rz(2.0430123) q[3];
sx q[3];
rz(-2.2618812) q[3];
sx q[3];
rz(-1.19899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
