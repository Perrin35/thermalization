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
rz(0.02778223) q[0];
sx q[0];
rz(1.8960928) q[0];
sx q[0];
rz(11.480042) q[0];
rz(-1.0983941) q[1];
sx q[1];
rz(-1.2135222) q[1];
sx q[1];
rz(-1.9953802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65685174) q[0];
sx q[0];
rz(-1.7324589) q[0];
sx q[0];
rz(-1.9009349) q[0];
x q[1];
rz(1.1623726) q[2];
sx q[2];
rz(-1.0193438) q[2];
sx q[2];
rz(-1.9806895) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0806345) q[1];
sx q[1];
rz(-0.74392156) q[1];
sx q[1];
rz(2.3638636) q[1];
rz(-pi) q[2];
rz(2.0729154) q[3];
sx q[3];
rz(-1.3716591) q[3];
sx q[3];
rz(0.56875436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(-2.7794465) q[2];
rz(-2.1892138) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(-2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1700965) q[0];
sx q[0];
rz(-2.9269452) q[0];
sx q[0];
rz(1.5317408) q[0];
rz(-1.8929298) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-0.85266399) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8448317) q[0];
sx q[0];
rz(-2.0728025) q[0];
sx q[0];
rz(0.39433483) q[0];
rz(-pi) q[1];
rz(1.8478099) q[2];
sx q[2];
rz(-1.4528034) q[2];
sx q[2];
rz(1.5867151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2679001) q[1];
sx q[1];
rz(-2.5136569) q[1];
sx q[1];
rz(1.2242641) q[1];
rz(2.144935) q[3];
sx q[3];
rz(-0.64145815) q[3];
sx q[3];
rz(1.5456089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93069211) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(-1.9453913) q[2];
rz(-3.1148552) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(-1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4836327) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(2.2657917) q[0];
rz(0.35975131) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(-2.13805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4205846) q[0];
sx q[0];
rz(-1.4996753) q[0];
sx q[0];
rz(2.5817602) q[0];
rz(-pi) q[1];
rz(-2.0802754) q[2];
sx q[2];
rz(-2.433314) q[2];
sx q[2];
rz(-2.2564611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0816312) q[1];
sx q[1];
rz(-1.7137032) q[1];
sx q[1];
rz(-2.4613454) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64522532) q[3];
sx q[3];
rz(-2.3564383) q[3];
sx q[3];
rz(-0.19715362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.935219) q[2];
sx q[2];
rz(-0.71722764) q[2];
sx q[2];
rz(2.3865872) q[2];
rz(-3.0435009) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(0.64495069) q[0];
rz(1.0954674) q[1];
sx q[1];
rz(-0.40830937) q[1];
sx q[1];
rz(-2.0300949) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43954014) q[0];
sx q[0];
rz(-1.5071535) q[0];
sx q[0];
rz(1.511277) q[0];
x q[1];
rz(-2.9745151) q[2];
sx q[2];
rz(-1.4896059) q[2];
sx q[2];
rz(2.2600391) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7202338) q[1];
sx q[1];
rz(-1.3610657) q[1];
sx q[1];
rz(-1.7214107) q[1];
rz(-pi) q[2];
rz(1.3502514) q[3];
sx q[3];
rz(-0.84201563) q[3];
sx q[3];
rz(0.33607863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5081386) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(0.48466361) q[2];
rz(-0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50705528) q[0];
sx q[0];
rz(-2.2978954) q[0];
sx q[0];
rz(1.7150568) q[0];
rz(2.2818395) q[1];
sx q[1];
rz(-2.8883002) q[1];
sx q[1];
rz(-2.3111129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62522307) q[0];
sx q[0];
rz(-0.88623669) q[0];
sx q[0];
rz(2.3806118) q[0];
x q[1];
rz(-0.057382345) q[2];
sx q[2];
rz(-2.2028366) q[2];
sx q[2];
rz(1.9789647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4791058) q[1];
sx q[1];
rz(-1.7104407) q[1];
sx q[1];
rz(-2.4100811) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4111341) q[3];
sx q[3];
rz(-1.9267907) q[3];
sx q[3];
rz(1.0245263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6972203) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(-0.74697906) q[2];
rz(0.30484453) q[3];
sx q[3];
rz(-1.7458785) q[3];
sx q[3];
rz(-1.9029118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759258) q[0];
sx q[0];
rz(-1.3015863) q[0];
sx q[0];
rz(-2.9081705) q[0];
rz(-0.4625136) q[1];
sx q[1];
rz(-1.1902483) q[1];
sx q[1];
rz(-0.40406427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44495917) q[0];
sx q[0];
rz(-0.6650266) q[0];
sx q[0];
rz(1.2103266) q[0];
rz(-pi) q[1];
rz(-2.7163804) q[2];
sx q[2];
rz(-1.8164486) q[2];
sx q[2];
rz(-0.52580357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3531472) q[1];
sx q[1];
rz(-2.569558) q[1];
sx q[1];
rz(-0.20010179) q[1];
rz(1.3498106) q[3];
sx q[3];
rz(-2.1897802) q[3];
sx q[3];
rz(-2.9429551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3290951) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(-2.4883032) q[2];
rz(1.4981883) q[3];
sx q[3];
rz(-0.61763063) q[3];
sx q[3];
rz(2.5352855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6752328) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(-0.32077041) q[0];
rz(2.4877211) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(3.0810862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9515647) q[0];
sx q[0];
rz(-1.4037644) q[0];
sx q[0];
rz(-1.8442979) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4662911) q[2];
sx q[2];
rz(-2.2847698) q[2];
sx q[2];
rz(1.7973822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.566129) q[1];
sx q[1];
rz(-0.82717973) q[1];
sx q[1];
rz(0.89218234) q[1];
rz(-pi) q[2];
rz(-2.7869446) q[3];
sx q[3];
rz(-1.186968) q[3];
sx q[3];
rz(0.5634748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(0.72236577) q[2];
rz(-2.8584621) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(-0.50962555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79001456) q[0];
sx q[0];
rz(-2.5836662) q[0];
sx q[0];
rz(0.18490069) q[0];
rz(0.50390538) q[1];
sx q[1];
rz(-1.3495219) q[1];
sx q[1];
rz(2.6769743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5123925) q[0];
sx q[0];
rz(-2.9058532) q[0];
sx q[0];
rz(2.3605973) q[0];
rz(-pi) q[1];
rz(-2.7242595) q[2];
sx q[2];
rz(-1.6586896) q[2];
sx q[2];
rz(-2.7869239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2459348) q[1];
sx q[1];
rz(-0.57809752) q[1];
sx q[1];
rz(-2.3996785) q[1];
rz(-pi) q[2];
rz(-2.9161386) q[3];
sx q[3];
rz(-2.4092393) q[3];
sx q[3];
rz(1.454753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23286143) q[2];
sx q[2];
rz(-1.5166914) q[2];
sx q[2];
rz(2.0880584) q[2];
rz(-2.4145224) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(-0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19186774) q[0];
sx q[0];
rz(-1.9504915) q[0];
sx q[0];
rz(-2.9560992) q[0];
rz(0.56974757) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(-1.1572908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13979736) q[0];
sx q[0];
rz(-1.9855864) q[0];
sx q[0];
rz(-1.7263361) q[0];
rz(-pi) q[1];
rz(0.050008372) q[2];
sx q[2];
rz(-0.93863025) q[2];
sx q[2];
rz(-0.061246733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2409199) q[1];
sx q[1];
rz(-1.0176588) q[1];
sx q[1];
rz(1.8227804) q[1];
rz(-pi) q[2];
rz(2.3367711) q[3];
sx q[3];
rz(-1.2535044) q[3];
sx q[3];
rz(0.99178335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40413228) q[2];
sx q[2];
rz(-0.22455939) q[2];
sx q[2];
rz(1.4004716) q[2];
rz(-2.776966) q[3];
sx q[3];
rz(-0.76321634) q[3];
sx q[3];
rz(-0.82272416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6309361) q[0];
sx q[0];
rz(-2.3031504) q[0];
sx q[0];
rz(3.1070218) q[0];
rz(-0.336054) q[1];
sx q[1];
rz(-2.2269109) q[1];
sx q[1];
rz(2.5539982) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97046554) q[0];
sx q[0];
rz(-1.1538527) q[0];
sx q[0];
rz(-1.4149026) q[0];
rz(0.35878344) q[2];
sx q[2];
rz(-2.3780895) q[2];
sx q[2];
rz(0.43649188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7149844) q[1];
sx q[1];
rz(-2.9054925) q[1];
sx q[1];
rz(1.3900422) q[1];
rz(-pi) q[2];
rz(-2.4617599) q[3];
sx q[3];
rz(-1.9380016) q[3];
sx q[3];
rz(2.1757371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33362886) q[2];
sx q[2];
rz(-2.2449292) q[2];
sx q[2];
rz(2.4133852) q[2];
rz(-2.7306469) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651466) q[0];
sx q[0];
rz(-1.5855753) q[0];
sx q[0];
rz(1.3517071) q[0];
rz(-2.7017055) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(2.7976703) q[2];
sx q[2];
rz(-1.2871598) q[2];
sx q[2];
rz(-3.1393307) q[2];
rz(-0.10671814) q[3];
sx q[3];
rz(-1.5244457) q[3];
sx q[3];
rz(2.5929858) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
