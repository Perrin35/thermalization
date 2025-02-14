OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(-1.902268) q[0];
sx q[0];
rz(-1.0898606) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(-0.90322948) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922674) q[0];
sx q[0];
rz(-0.9454596) q[0];
sx q[0];
rz(2.3614285) q[0];
rz(0.014711424) q[2];
sx q[2];
rz(-0.88764578) q[2];
sx q[2];
rz(-2.8840182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0923903) q[1];
sx q[1];
rz(-0.70194879) q[1];
sx q[1];
rz(-0.54767056) q[1];
rz(-1.0711477) q[3];
sx q[3];
rz(-1.6766492) q[3];
sx q[3];
rz(-1.203804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.696306) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(2.7083) q[2];
rz(-2.8516234) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(2.8210631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1746154) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(1.2906661) q[0];
rz(-0.67944747) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(-2.2490833) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9226869) q[0];
sx q[0];
rz(-1.3935879) q[0];
sx q[0];
rz(-0.025434504) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65373924) q[2];
sx q[2];
rz(-1.8934665) q[2];
sx q[2];
rz(-1.309909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1074236) q[1];
sx q[1];
rz(-2.4236107) q[1];
sx q[1];
rz(-0.42426829) q[1];
rz(-pi) q[2];
rz(-2.1492747) q[3];
sx q[3];
rz(-2.7035993) q[3];
sx q[3];
rz(-2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14900011) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(-1.1326257) q[2];
rz(-2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76729524) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(2.1597916) q[0];
rz(-0.94217316) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(2.2312677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8273329) q[0];
sx q[0];
rz(-1.9935441) q[0];
sx q[0];
rz(-2.6092922) q[0];
x q[1];
rz(-0.40159638) q[2];
sx q[2];
rz(-2.842428) q[2];
sx q[2];
rz(-2.6515692) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83214789) q[1];
sx q[1];
rz(-2.6845686) q[1];
sx q[1];
rz(2.3085269) q[1];
x q[2];
rz(1.0873454) q[3];
sx q[3];
rz(-0.48763613) q[3];
sx q[3];
rz(-0.82043649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(-2.9777891) q[3];
sx q[3];
rz(-0.87023321) q[3];
sx q[3];
rz(0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(-0.69766587) q[0];
rz(0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(1.9465416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9126667) q[0];
sx q[0];
rz(-0.9879092) q[0];
sx q[0];
rz(2.7253662) q[0];
rz(-pi) q[1];
rz(0.61364737) q[2];
sx q[2];
rz(-0.23229182) q[2];
sx q[2];
rz(0.42343806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6704441) q[1];
sx q[1];
rz(-0.19757825) q[1];
sx q[1];
rz(-2.1414645) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8175587) q[3];
sx q[3];
rz(-0.29280782) q[3];
sx q[3];
rz(-0.68835252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-2.9441693) q[2];
rz(0.10489634) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(-3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5474434) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(-1.0092258) q[0];
rz(-3.1127473) q[1];
sx q[1];
rz(-1.6639158) q[1];
sx q[1];
rz(2.4829594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65153367) q[0];
sx q[0];
rz(-2.5680827) q[0];
sx q[0];
rz(-2.308368) q[0];
rz(1.0844857) q[2];
sx q[2];
rz(-1.3895814) q[2];
sx q[2];
rz(1.043037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11396961) q[1];
sx q[1];
rz(-1.4024892) q[1];
sx q[1];
rz(1.6644003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4464408) q[3];
sx q[3];
rz(-1.3926695) q[3];
sx q[3];
rz(-3.0722838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7625526) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(-0.23400447) q[2];
rz(1.0903357) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(-1.0696629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91319084) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(1.2983904) q[0];
rz(2.3550854) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(1.7826805) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494736) q[0];
sx q[0];
rz(-0.91128659) q[0];
sx q[0];
rz(-1.6716206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26769079) q[2];
sx q[2];
rz(-1.7867733) q[2];
sx q[2];
rz(-2.9962073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.003215) q[1];
sx q[1];
rz(-1.5724465) q[1];
sx q[1];
rz(-1.5720075) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5868191) q[3];
sx q[3];
rz(-2.6027461) q[3];
sx q[3];
rz(-0.19308819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97335029) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(2.8016222) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.2350157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17609734) q[0];
sx q[0];
rz(-1.9416932) q[0];
sx q[0];
rz(-2.4489016) q[0];
x q[1];
rz(2.6831362) q[2];
sx q[2];
rz(-1.2453015) q[2];
sx q[2];
rz(2.1874962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75720471) q[1];
sx q[1];
rz(-1.5761307) q[1];
sx q[1];
rz(2.6621006) q[1];
rz(-pi) q[2];
x q[2];
rz(1.791802) q[3];
sx q[3];
rz(-1.955893) q[3];
sx q[3];
rz(2.4996595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0299224) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(2.5862582) q[2];
rz(-0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.70917201) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-1.1093371) q[0];
rz(-2.0093911) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(-2.1999377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718194) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(0.18751796) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(2.6251453) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.718457) q[1];
sx q[1];
rz(-0.65538156) q[1];
sx q[1];
rz(1.7792367) q[1];
rz(1.7835983) q[3];
sx q[3];
rz(-1.1929242) q[3];
sx q[3];
rz(1.4546284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7534916) q[2];
sx q[2];
rz(-2.1145623) q[2];
sx q[2];
rz(1.5388185) q[2];
rz(1.3711035) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(-3.0901618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082315363) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(2.1737461) q[0];
rz(1.8702501) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4000452) q[0];
sx q[0];
rz(-1.1843268) q[0];
sx q[0];
rz(2.5611112) q[0];
rz(-pi) q[1];
rz(-0.25904663) q[2];
sx q[2];
rz(-1.7648089) q[2];
sx q[2];
rz(1.8679179) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4092952) q[1];
sx q[1];
rz(-1.4753072) q[1];
sx q[1];
rz(0.68640253) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8621919) q[3];
sx q[3];
rz(-1.8899922) q[3];
sx q[3];
rz(0.58871692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1900078) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(1.9742981) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(-2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9645914) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(-0.28468537) q[0];
rz(3.0668861) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(1.7237192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95686326) q[0];
sx q[0];
rz(-1.7635582) q[0];
sx q[0];
rz(-2.3817029) q[0];
rz(-1.013834) q[2];
sx q[2];
rz(-1.6979453) q[2];
sx q[2];
rz(-0.47455041) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26894179) q[1];
sx q[1];
rz(-1.1564019) q[1];
sx q[1];
rz(0.68434836) q[1];
rz(-pi) q[2];
rz(-1.5261493) q[3];
sx q[3];
rz(-0.64638019) q[3];
sx q[3];
rz(-0.97585362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(-1.2791862) q[2];
rz(2.159436) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(2.2568683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2411156) q[0];
sx q[0];
rz(-1.3237088) q[0];
sx q[0];
rz(1.5644912) q[0];
rz(-1.0152394) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(-1.3080636) q[2];
sx q[2];
rz(-1.4663525) q[2];
sx q[2];
rz(2.732085) q[2];
rz(0.61839907) q[3];
sx q[3];
rz(-2.3000345) q[3];
sx q[3];
rz(-2.2207501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
