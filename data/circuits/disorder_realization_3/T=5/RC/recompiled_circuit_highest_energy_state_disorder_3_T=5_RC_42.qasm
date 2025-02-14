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
rz(-1.2454998) q[0];
sx q[0];
rz(1.0863289) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(-1.9280704) q[1];
sx q[1];
rz(1.9953802) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47482309) q[0];
sx q[0];
rz(-2.7753029) q[0];
sx q[0];
rz(1.1046871) q[0];
rz(-pi) q[1];
rz(-0.59046794) q[2];
sx q[2];
rz(-1.2256978) q[2];
sx q[2];
rz(-2.9546628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0806345) q[1];
sx q[1];
rz(-2.3976711) q[1];
sx q[1];
rz(2.3638636) q[1];
rz(-pi) q[2];
rz(-2.0729154) q[3];
sx q[3];
rz(-1.7699336) q[3];
sx q[3];
rz(0.56875436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9132797) q[2];
sx q[2];
rz(-1.5659119) q[2];
sx q[2];
rz(-2.7794465) q[2];
rz(-2.1892138) q[3];
sx q[3];
rz(-2.5824472) q[3];
sx q[3];
rz(2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1700965) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(1.6098518) q[0];
rz(1.8929298) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-2.2889287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.296761) q[0];
sx q[0];
rz(-2.0728025) q[0];
sx q[0];
rz(-2.7472578) q[0];
x q[1];
rz(0.12262362) q[2];
sx q[2];
rz(-1.8458335) q[2];
sx q[2];
rz(-0.017539121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2679001) q[1];
sx q[1];
rz(-2.5136569) q[1];
sx q[1];
rz(-1.9173286) q[1];
rz(2.144935) q[3];
sx q[3];
rz(-0.64145815) q[3];
sx q[3];
rz(-1.5959838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2109005) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(1.1962013) q[2];
rz(-3.1148552) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(1.4727717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4836327) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(-0.87580097) q[0];
rz(0.35975131) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(1.0035427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7210081) q[0];
sx q[0];
rz(-1.4996753) q[0];
sx q[0];
rz(2.5817602) q[0];
rz(-pi) q[1];
rz(0.92873145) q[2];
sx q[2];
rz(-1.8936529) q[2];
sx q[2];
rz(-1.0869458) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4565125) q[1];
sx q[1];
rz(-2.4488423) q[1];
sx q[1];
rz(-2.9167006) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-2.935219) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(-2.3865872) q[2];
rz(0.098091789) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(-0.64495069) q[0];
rz(-1.0954674) q[1];
sx q[1];
rz(-0.40830937) q[1];
sx q[1];
rz(-1.1114978) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0065466) q[0];
sx q[0];
rz(-1.5113976) q[0];
sx q[0];
rz(-3.0778372) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16707755) q[2];
sx q[2];
rz(-1.4896059) q[2];
sx q[2];
rz(-0.88155356) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9330628) q[1];
sx q[1];
rz(-0.25756076) q[1];
sx q[1];
rz(2.527586) q[1];
rz(-pi) q[2];
rz(-0.74097775) q[3];
sx q[3];
rz(-1.40687) q[3];
sx q[3];
rz(1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5081386) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(-0.48466361) q[2];
rz(0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6345374) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(1.7150568) q[0];
rz(0.85975319) q[1];
sx q[1];
rz(-2.8883002) q[1];
sx q[1];
rz(2.3111129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4035506) q[0];
sx q[0];
rz(-2.134424) q[0];
sx q[0];
rz(-0.72569816) q[0];
rz(-3.0842103) q[2];
sx q[2];
rz(-2.2028366) q[2];
sx q[2];
rz(1.1626279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0326091) q[1];
sx q[1];
rz(-2.2936037) q[1];
sx q[1];
rz(-1.384114) q[1];
rz(0.36019616) q[3];
sx q[3];
rz(-1.4212228) q[3];
sx q[3];
rz(0.60233483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44437235) q[2];
sx q[2];
rz(-2.1264919) q[2];
sx q[2];
rz(-0.74697906) q[2];
rz(2.8367481) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9656669) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(2.9081705) q[0];
rz(-2.6790791) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(2.7375284) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44495917) q[0];
sx q[0];
rz(-2.4765661) q[0];
sx q[0];
rz(-1.9312661) q[0];
rz(-pi) q[1];
rz(0.42521221) q[2];
sx q[2];
rz(-1.325144) q[2];
sx q[2];
rz(0.52580357) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.048745884) q[1];
sx q[1];
rz(-1.6786075) q[1];
sx q[1];
rz(-2.5786933) q[1];
x q[2];
rz(1.791782) q[3];
sx q[3];
rz(-0.95181247) q[3];
sx q[3];
rz(0.19863752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3290951) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(-2.4883032) q[2];
rz(-1.6434044) q[3];
sx q[3];
rz(-2.523962) q[3];
sx q[3];
rz(-2.5352855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752328) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(0.32077041) q[0];
rz(2.4877211) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(3.0810862) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.190028) q[0];
sx q[0];
rz(-1.7378283) q[0];
sx q[0];
rz(1.2972948) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73326335) q[2];
sx q[2];
rz(-1.0787233) q[2];
sx q[2];
rz(0.70962188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.566129) q[1];
sx q[1];
rz(-2.3144129) q[1];
sx q[1];
rz(-0.89218234) q[1];
x q[2];
rz(-0.860608) q[3];
sx q[3];
rz(-0.51653701) q[3];
sx q[3];
rz(-2.9252657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(-2.4192269) q[2];
rz(-0.28313053) q[3];
sx q[3];
rz(-0.88909283) q[3];
sx q[3];
rz(-0.50962555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515781) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(-0.18490069) q[0];
rz(2.6376873) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(-0.46461836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9757742) q[0];
sx q[0];
rz(-1.7374455) q[0];
sx q[0];
rz(-1.7383132) q[0];
rz(-pi) q[1];
rz(-1.6668929) q[2];
sx q[2];
rz(-1.155174) q[2];
sx q[2];
rz(1.1772275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89565784) q[1];
sx q[1];
rz(-0.57809752) q[1];
sx q[1];
rz(-0.74191414) q[1];
rz(-0.71962756) q[3];
sx q[3];
rz(-1.4207645) q[3];
sx q[3];
rz(0.052879083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23286143) q[2];
sx q[2];
rz(-1.5166914) q[2];
sx q[2];
rz(-2.0880584) q[2];
rz(0.72707027) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(2.5963636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19186774) q[0];
sx q[0];
rz(-1.9504915) q[0];
sx q[0];
rz(-2.9560992) q[0];
rz(2.5718451) q[1];
sx q[1];
rz(-2.2173917) q[1];
sx q[1];
rz(-1.1572908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3678903) q[0];
sx q[0];
rz(-1.4285402) q[0];
sx q[0];
rz(-0.41928798) q[0];
rz(-0.93803381) q[2];
sx q[2];
rz(-1.530458) q[2];
sx q[2];
rz(-1.4799839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9461575) q[1];
sx q[1];
rz(-1.7845673) q[1];
sx q[1];
rz(2.5740088) q[1];
x q[2];
rz(-0.42753856) q[3];
sx q[3];
rz(-0.85179177) q[3];
sx q[3];
rz(2.2710272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40413228) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(-1.4004716) q[2];
rz(0.36462668) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(0.82272416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.51065651) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(-0.034570845) q[0];
rz(-0.336054) q[1];
sx q[1];
rz(-0.91468179) q[1];
sx q[1];
rz(-2.5539982) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776992) q[0];
sx q[0];
rz(-1.4283533) q[0];
sx q[0];
rz(0.42147984) q[0];
x q[1];
rz(-2.7828092) q[2];
sx q[2];
rz(-0.76350313) q[2];
sx q[2];
rz(-0.43649188) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61240367) q[1];
sx q[1];
rz(-1.3386139) q[1];
sx q[1];
rz(0.043223765) q[1];
rz(-pi) q[2];
rz(0.54909191) q[3];
sx q[3];
rz(-2.3830416) q[3];
sx q[3];
rz(-2.9545934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8079638) q[2];
sx q[2];
rz(-2.2449292) q[2];
sx q[2];
rz(-2.4133852) q[2];
rz(0.41094574) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644607) q[0];
sx q[0];
rz(-1.5560173) q[0];
sx q[0];
rz(-1.7898855) q[0];
rz(0.43988718) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(-2.7976703) q[2];
sx q[2];
rz(-1.8544329) q[2];
sx q[2];
rz(0.0022619958) q[2];
rz(1.6174118) q[3];
sx q[3];
rz(-1.4641932) q[3];
sx q[3];
rz(-2.1144397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
