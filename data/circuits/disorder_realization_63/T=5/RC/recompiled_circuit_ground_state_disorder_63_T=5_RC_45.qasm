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
rz(1.0898606) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(-1.8884594) q[1];
sx q[1];
rz(0.90322948) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3147557) q[0];
sx q[0];
rz(-2.1776548) q[0];
sx q[0];
rz(0.77758247) q[0];
rz(3.1268812) q[2];
sx q[2];
rz(-0.88764578) q[2];
sx q[2];
rz(2.8840182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62476872) q[1];
sx q[1];
rz(-0.98691578) q[1];
sx q[1];
rz(1.9855687) q[1];
x q[2];
rz(2.0704449) q[3];
sx q[3];
rz(-1.4649434) q[3];
sx q[3];
rz(-1.9377886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.696306) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(0.43329263) q[2];
rz(2.8516234) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(-0.32052952) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9669773) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(1.2906661) q[0];
rz(0.67944747) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(0.89250934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941867) q[0];
sx q[0];
rz(-1.5958324) q[0];
sx q[0];
rz(1.3935318) q[0];
rz(2.638916) q[2];
sx q[2];
rz(-2.423175) q[2];
sx q[2];
rz(0.65332149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50607562) q[1];
sx q[1];
rz(-2.2137224) q[1];
sx q[1];
rz(1.2256114) q[1];
x q[2];
rz(0.99231798) q[3];
sx q[3];
rz(-2.7035993) q[3];
sx q[3];
rz(0.43644825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9925925) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(-2.0089669) q[2];
rz(0.66631404) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(-1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3742974) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(2.1597916) q[0];
rz(-0.94217316) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(2.2312677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019429723) q[0];
sx q[0];
rz(-1.0896026) q[0];
sx q[0];
rz(-1.0896171) q[0];
rz(-pi) q[1];
rz(1.6907755) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-0.071823013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.040346216) q[1];
sx q[1];
rz(-1.2381499) q[1];
sx q[1];
rz(0.31942792) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0542473) q[3];
sx q[3];
rz(-0.48763613) q[3];
sx q[3];
rz(-2.3211562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(0.18690898) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-0.87023321) q[3];
sx q[3];
rz(-0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1214889) q[0];
sx q[0];
rz(-1.4153471) q[0];
sx q[0];
rz(-2.4439268) q[0];
rz(0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.9465416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89681399) q[0];
sx q[0];
rz(-1.2264826) q[0];
sx q[0];
rz(-2.1953775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19104345) q[2];
sx q[2];
rz(-1.7037539) q[2];
sx q[2];
rz(1.7482479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89151284) q[1];
sx q[1];
rz(-1.7367474) q[1];
sx q[1];
rz(-0.1077229) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2863808) q[3];
sx q[3];
rz(-1.5002325) q[3];
sx q[3];
rz(-2.0224935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-0.19742337) q[2];
rz(0.10489634) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(0.088002861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65153367) q[0];
sx q[0];
rz(-2.5680827) q[0];
sx q[0];
rz(2.308368) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20435996) q[2];
sx q[2];
rz(-2.0484701) q[2];
sx q[2];
rz(0.43276873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4411021) q[1];
sx q[1];
rz(-1.6630739) q[1];
sx q[1];
rz(2.9725595) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8010718) q[3];
sx q[3];
rz(-0.88876969) q[3];
sx q[3];
rz(1.7868228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7625526) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(-0.23400447) q[2];
rz(-2.0512569) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(2.0719297) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2284018) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(-1.3589121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14058622) q[0];
sx q[0];
rz(-1.4911665) q[0];
sx q[0];
rz(-2.4796159) q[0];
rz(-pi) q[1];
rz(2.8739019) q[2];
sx q[2];
rz(-1.7867733) q[2];
sx q[2];
rz(-2.9962073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1383776) q[1];
sx q[1];
rz(-1.5691461) q[1];
sx q[1];
rz(-1.5720075) q[1];
rz(0.5868191) q[3];
sx q[3];
rz(-0.53884655) q[3];
sx q[3];
rz(-2.9485045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(-2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(-0.33997047) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(-1.9065769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.039073) q[0];
sx q[0];
rz(-2.2082941) q[0];
sx q[0];
rz(2.0387285) q[0];
rz(-1.9307617) q[2];
sx q[2];
rz(-1.1380929) q[2];
sx q[2];
rz(-0.77322799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81081796) q[1];
sx q[1];
rz(-2.0502809) q[1];
sx q[1];
rz(-1.5768087) q[1];
rz(-1.3497906) q[3];
sx q[3];
rz(-1.1856996) q[3];
sx q[3];
rz(-2.4996595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(2.5862582) q[2];
rz(0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(-0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324206) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-1.1093371) q[0];
rz(2.0093911) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(0.94165492) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8245715) q[0];
sx q[0];
rz(-0.57972687) q[0];
sx q[0];
rz(1.8648022) q[0];
rz(-pi) q[1];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(2.6251453) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98151614) q[1];
sx q[1];
rz(-1.6972516) q[1];
sx q[1];
rz(-0.92595788) q[1];
rz(-1.7835983) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(-1.6869643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(-1.5388185) q[2];
rz(1.3711035) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082315363) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(0.96784651) q[0];
rz(-1.2713426) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74154749) q[0];
sx q[0];
rz(-1.1843268) q[0];
sx q[0];
rz(-2.5611112) q[0];
rz(-pi) q[1];
rz(-2.4872753) q[2];
sx q[2];
rz(-2.8192602) q[2];
sx q[2];
rz(2.2155025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3810515) q[1];
sx q[1];
rz(-0.88812056) q[1];
sx q[1];
rz(-1.447601) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8621919) q[3];
sx q[3];
rz(-1.8899922) q[3];
sx q[3];
rz(-2.5528757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(-1.1672945) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.17700125) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(-2.8569073) q[0];
rz(-0.074706569) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(-1.4178735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79399949) q[0];
sx q[0];
rz(-0.82836223) q[0];
sx q[0];
rz(-1.8338127) q[0];
rz(-pi) q[1];
rz(-1.333513) q[2];
sx q[2];
rz(-2.5717989) q[2];
sx q[2];
rz(-2.2462318) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3813044) q[1];
sx q[1];
rz(-2.3592296) q[1];
sx q[1];
rz(2.5336877) q[1];
x q[2];
rz(-3.10793) q[3];
sx q[3];
rz(-2.2164248) q[3];
sx q[3];
rz(2.1098304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(1.2791862) q[2];
rz(-0.98215669) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(0.88472432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2411156) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(1.0152394) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(1.1871836) q[2];
sx q[2];
rz(-2.8593079) q[2];
sx q[2];
rz(1.5310892) q[2];
rz(-2.5231936) q[3];
sx q[3];
rz(-2.3000345) q[3];
sx q[3];
rz(-2.2207501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
