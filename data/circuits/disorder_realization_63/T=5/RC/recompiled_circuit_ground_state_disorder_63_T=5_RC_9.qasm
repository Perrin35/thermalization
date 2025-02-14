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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82683691) q[0];
sx q[0];
rz(-0.96393782) q[0];
sx q[0];
rz(-2.3640102) q[0];
rz(-pi) q[1];
rz(0.8875928) q[2];
sx q[2];
rz(-1.5593865) q[2];
sx q[2];
rz(-1.3225088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1841539) q[1];
sx q[1];
rz(-1.9136962) q[1];
sx q[1];
rz(-2.5162906) q[1];
rz(-pi) q[2];
rz(3.0211307) q[3];
sx q[3];
rz(-1.0742004) q[3];
sx q[3];
rz(2.7169926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4452867) q[2];
sx q[2];
rz(-0.13662766) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(-0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-2.8210631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1746154) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(-1.2906661) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(-0.89250934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941867) q[0];
sx q[0];
rz(-1.5457602) q[0];
sx q[0];
rz(1.7480609) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50267668) q[2];
sx q[2];
rz(-2.423175) q[2];
sx q[2];
rz(2.4882712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50607562) q[1];
sx q[1];
rz(-2.2137224) q[1];
sx q[1];
rz(-1.2256114) q[1];
rz(-pi) q[2];
rz(1.1970911) q[3];
sx q[3];
rz(-1.804816) q[3];
sx q[3];
rz(-1.6683874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14900011) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(-1.1326257) q[2];
rz(-2.4752786) q[3];
sx q[3];
rz(-1.7814813) q[3];
sx q[3];
rz(-1.2247491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.1994195) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(-0.91032496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1221629) q[0];
sx q[0];
rz(-1.0896026) q[0];
sx q[0];
rz(-1.0896171) q[0];
rz(-pi) q[1];
rz(0.40159638) q[2];
sx q[2];
rz(-0.29916468) q[2];
sx q[2];
rz(-2.6515692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7187319) q[1];
sx q[1];
rz(-1.2694468) q[1];
sx q[1];
rz(1.919793) q[1];
x q[2];
rz(0.24170928) q[3];
sx q[3];
rz(-1.998566) q[3];
sx q[3];
rz(0.28423968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(-0.16380353) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201037) q[0];
sx q[0];
rz(-1.4153471) q[0];
sx q[0];
rz(-0.69766587) q[0];
rz(-0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(-1.9465416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9055525) q[0];
sx q[0];
rz(-0.70193203) q[0];
sx q[0];
rz(1.0206971) q[0];
rz(-0.19104345) q[2];
sx q[2];
rz(-1.7037539) q[2];
sx q[2];
rz(-1.7482479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4711485) q[1];
sx q[1];
rz(-2.9440144) q[1];
sx q[1];
rz(-2.1414645) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8552119) q[3];
sx q[3];
rz(-1.6413601) q[3];
sx q[3];
rz(-1.1190991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4635072) q[2];
sx q[2];
rz(-1.4231851) q[2];
sx q[2];
rz(0.19742337) q[2];
rz(3.0366963) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5941493) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(2.1323668) q[0];
rz(3.1127473) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(-0.65863329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26732609) q[0];
sx q[0];
rz(-1.1972885) q[0];
sx q[0];
rz(-1.1248571) q[0];
rz(-1.1971742) q[2];
sx q[2];
rz(-2.6251617) q[2];
sx q[2];
rz(2.2852798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62477905) q[1];
sx q[1];
rz(-2.9492231) q[1];
sx q[1];
rz(-2.6386847) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69515184) q[3];
sx q[3];
rz(-1.7489232) q[3];
sx q[3];
rz(3.0722838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7625526) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(2.9075882) q[2];
rz(-1.0903357) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91319084) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(1.7826805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8131065) q[0];
sx q[0];
rz(-0.66603249) q[0];
sx q[0];
rz(-3.0124927) q[0];
rz(-0.26769079) q[2];
sx q[2];
rz(-1.3548193) q[2];
sx q[2];
rz(2.9962073) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56758331) q[1];
sx q[1];
rz(-1.5695851) q[1];
sx q[1];
rz(3.1399425) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6796713) q[3];
sx q[3];
rz(-1.2826903) q[3];
sx q[3];
rz(1.896331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4141554) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(-2.9690913) q[2];
rz(-0.48834673) q[3];
sx q[3];
rz(-1.9988873) q[3];
sx q[3];
rz(-1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97335029) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(-0.33997047) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.2350157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8064855) q[0];
sx q[0];
rz(-0.77101427) q[0];
sx q[0];
rz(2.5946027) q[0];
x q[1];
rz(-1.9307617) q[2];
sx q[2];
rz(-2.0034997) q[2];
sx q[2];
rz(-2.3683647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3177424) q[1];
sx q[1];
rz(-2.6620733) q[1];
sx q[1];
rz(3.13003) q[1];
rz(0.49576393) q[3];
sx q[3];
rz(-0.44124441) q[3];
sx q[3];
rz(3.038681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(2.9739144) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(-2.663747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324206) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(2.0322556) q[0];
rz(-2.0093911) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(0.94165492) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718194) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(-2.9540747) q[0];
x q[1];
rz(2.9865725) q[2];
sx q[2];
rz(-1.2885254) q[2];
sx q[2];
rz(2.043743) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.718457) q[1];
sx q[1];
rz(-2.4862111) q[1];
sx q[1];
rz(1.7792367) q[1];
rz(-pi) q[2];
rz(-1.3579943) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(1.6869643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38810101) q[2];
sx q[2];
rz(-2.1145623) q[2];
sx q[2];
rz(-1.6027742) q[2];
rz(1.7704891) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(3.0901618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082315363) q[0];
sx q[0];
rz(-2.5264854) q[0];
sx q[0];
rz(-0.96784651) q[0];
rz(1.8702501) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(-0.06180067) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58690155) q[0];
sx q[0];
rz(-2.1036316) q[0];
sx q[0];
rz(-2.0237049) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.882546) q[2];
sx q[2];
rz(-1.3767837) q[2];
sx q[2];
rz(-1.2736748) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1872431) q[1];
sx q[1];
rz(-0.69194416) q[1];
sx q[1];
rz(-0.14999565) q[1];
rz(-pi) q[2];
rz(-1.8621919) q[3];
sx q[3];
rz(-1.2516004) q[3];
sx q[3];
rz(-0.58871692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9515848) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(-1.9742981) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41499785) q[0];
sx q[0];
rz(-2.3624067) q[0];
sx q[0];
rz(-0.27611538) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1277587) q[2];
sx q[2];
rz(-1.4436473) q[2];
sx q[2];
rz(0.47455041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7602882) q[1];
sx q[1];
rz(-0.78236303) q[1];
sx q[1];
rz(-2.5336877) q[1];
x q[2];
rz(2.2166972) q[3];
sx q[3];
rz(-1.5439111) q[3];
sx q[3];
rz(0.55929375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3384) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(-1.8624064) q[2];
rz(2.159436) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(2.2568683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(2.1263532) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(-1.954409) q[2];
sx q[2];
rz(-2.8593079) q[2];
sx q[2];
rz(1.5310892) q[2];
rz(0.99526631) q[3];
sx q[3];
rz(-2.223816) q[3];
sx q[3];
rz(1.73903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
