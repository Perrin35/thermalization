OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(-2.0517321) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(2.2383632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21891864) q[0];
sx q[0];
rz(-2.196133) q[0];
sx q[0];
rz(2.3614285) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2539999) q[2];
sx q[2];
rz(-1.5822062) q[2];
sx q[2];
rz(1.8190839) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62476872) q[1];
sx q[1];
rz(-0.98691578) q[1];
sx q[1];
rz(1.1560239) q[1];
rz(-pi) q[2];
rz(1.3525659) q[3];
sx q[3];
rz(-0.50980836) q[3];
sx q[3];
rz(-2.9657983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(-2.8516234) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9669773) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(-1.8509266) q[0];
rz(-2.4621452) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(2.2490833) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21890579) q[0];
sx q[0];
rz(-1.7480047) q[0];
sx q[0];
rz(-0.025434504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1721482) q[2];
sx q[2];
rz(-0.95602334) q[2];
sx q[2];
rz(3.1190256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.034169056) q[1];
sx q[1];
rz(-0.71798199) q[1];
sx q[1];
rz(-2.7173244) q[1];
x q[2];
rz(2.1492747) q[3];
sx q[3];
rz(-0.43799339) q[3];
sx q[3];
rz(0.43644825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9925925) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(2.0089669) q[2];
rz(-2.4752786) q[3];
sx q[3];
rz(-1.7814813) q[3];
sx q[3];
rz(1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742974) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(0.98180109) q[0];
rz(-2.1994195) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(2.2312677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8273329) q[0];
sx q[0];
rz(-1.9935441) q[0];
sx q[0];
rz(2.6092922) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40159638) q[2];
sx q[2];
rz(-2.842428) q[2];
sx q[2];
rz(-0.49002346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83214789) q[1];
sx q[1];
rz(-2.6845686) q[1];
sx q[1];
rz(0.8330658) q[1];
x q[2];
rz(0.24170928) q[3];
sx q[3];
rz(-1.998566) q[3];
sx q[3];
rz(0.28423968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(0.18690898) q[2];
rz(-0.16380353) q[3];
sx q[3];
rz(-0.87023321) q[3];
sx q[3];
rz(3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214889) q[0];
sx q[0];
rz(-1.4153471) q[0];
sx q[0];
rz(-0.69766587) q[0];
rz(2.5949219) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.1950511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.228926) q[0];
sx q[0];
rz(-2.1536835) q[0];
sx q[0];
rz(-0.41622644) q[0];
x q[1];
rz(2.9505492) q[2];
sx q[2];
rz(-1.4378387) q[2];
sx q[2];
rz(-1.3933448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6704441) q[1];
sx q[1];
rz(-2.9440144) q[1];
sx q[1];
rz(1.0001282) q[1];
rz(-pi) q[2];
rz(3.0680858) q[3];
sx q[3];
rz(-1.8544844) q[3];
sx q[3];
rz(2.7105041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-2.9441693) q[2];
rz(-3.0366963) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(0.088002861) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(2.1323668) q[0];
rz(3.1127473) q[1];
sx q[1];
rz(-1.6639158) q[1];
sx q[1];
rz(-2.4829594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26732609) q[0];
sx q[0];
rz(-1.1972885) q[0];
sx q[0];
rz(1.1248571) q[0];
x q[1];
rz(1.9444185) q[2];
sx q[2];
rz(-2.6251617) q[2];
sx q[2];
rz(2.2852798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4411021) q[1];
sx q[1];
rz(-1.6630739) q[1];
sx q[1];
rz(0.16903318) q[1];
x q[2];
rz(0.27401383) q[3];
sx q[3];
rz(-2.4276795) q[3];
sx q[3];
rz(1.4307724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.37904) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(2.9075882) q[2];
rz(-2.0512569) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(-2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2284018) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(-1.2983904) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(1.3589121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0010064) q[0];
sx q[0];
rz(-1.4911665) q[0];
sx q[0];
rz(0.66197673) q[0];
rz(-pi) q[1];
rz(-1.7944906) q[2];
sx q[2];
rz(-1.8321206) q[2];
sx q[2];
rz(1.4841207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.003215) q[1];
sx q[1];
rz(-1.5691461) q[1];
sx q[1];
rz(-1.5720075) q[1];
rz(-pi) q[2];
rz(0.5868191) q[3];
sx q[3];
rz(-0.53884655) q[3];
sx q[3];
rz(0.19308819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7274373) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(-2.9690913) q[2];
rz(0.48834673) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-1.1138227) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1682424) q[0];
sx q[0];
rz(-2.2689447) q[0];
sx q[0];
rz(2.8016222) q[0];
rz(-2.8151457) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(1.2350157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9654953) q[0];
sx q[0];
rz(-1.1998994) q[0];
sx q[0];
rz(0.69269106) q[0];
rz(-2.6831362) q[2];
sx q[2];
rz(-1.2453015) q[2];
sx q[2];
rz(0.95409648) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82385027) q[1];
sx q[1];
rz(-2.6620733) q[1];
sx q[1];
rz(-3.13003) q[1];
rz(-1.3497906) q[3];
sx q[3];
rz(-1.1856996) q[3];
sx q[3];
rz(0.64193314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0299224) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(2.663747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324206) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(2.0093911) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(2.1999377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6397259) q[0];
sx q[0];
rz(-1.7302156) q[0];
sx q[0];
rz(1.0109883) q[0];
rz(-pi) q[1];
rz(-0.15502013) q[2];
sx q[2];
rz(-1.2885254) q[2];
sx q[2];
rz(-1.0978497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1600765) q[1];
sx q[1];
rz(-1.4443411) q[1];
sx q[1];
rz(0.92595788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48896472) q[3];
sx q[3];
rz(-0.43114907) q[3];
sx q[3];
rz(-2.2167689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7534916) q[2];
sx q[2];
rz(-2.1145623) q[2];
sx q[2];
rz(-1.6027742) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(3.0901618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(-2.1737461) q[0];
rz(1.8702501) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4000452) q[0];
sx q[0];
rz(-1.9572659) q[0];
sx q[0];
rz(-0.58048141) q[0];
rz(-pi) q[1];
x q[1];
rz(2.882546) q[2];
sx q[2];
rz(-1.7648089) q[2];
sx q[2];
rz(-1.2736748) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73229746) q[1];
sx q[1];
rz(-1.4753072) q[1];
sx q[1];
rz(-0.68640253) q[1];
x q[2];
rz(-1.2794007) q[3];
sx q[3];
rz(-1.2516004) q[3];
sx q[3];
rz(0.58871692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(-3.0657213) q[2];
rz(-1.1672945) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.17700125) q[0];
sx q[0];
rz(-0.27934203) q[0];
sx q[0];
rz(-0.28468537) q[0];
rz(0.074706569) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(1.4178735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41499785) q[0];
sx q[0];
rz(-0.779186) q[0];
sx q[0];
rz(0.27611538) q[0];
x q[1];
rz(2.1277587) q[2];
sx q[2];
rz(-1.6979453) q[2];
sx q[2];
rz(2.6670422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26894179) q[1];
sx q[1];
rz(-1.9851908) q[1];
sx q[1];
rz(-0.68434836) q[1];
x q[2];
rz(0.92489544) q[3];
sx q[3];
rz(-1.5439111) q[3];
sx q[3];
rz(2.5822989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3384) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(-1.8624064) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(-2.2568683) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(-2.1263532) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(1.8335291) q[2];
sx q[2];
rz(-1.4663525) q[2];
sx q[2];
rz(2.732085) q[2];
rz(-0.99526631) q[3];
sx q[3];
rz(-0.91777663) q[3];
sx q[3];
rz(-1.4025626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
