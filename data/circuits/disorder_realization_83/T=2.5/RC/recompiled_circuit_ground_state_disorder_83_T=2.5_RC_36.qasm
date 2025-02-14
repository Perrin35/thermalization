OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72252005) q[0];
sx q[0];
rz(-0.68113405) q[0];
sx q[0];
rz(-1.0406159) q[0];
rz(2.430727) q[1];
sx q[1];
rz(3.791888) q[1];
sx q[1];
rz(11.314582) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1052763) q[0];
sx q[0];
rz(-1.9364089) q[0];
sx q[0];
rz(2.7514639) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5898143) q[2];
sx q[2];
rz(-1.1349196) q[2];
sx q[2];
rz(2.9575155) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2179778) q[1];
sx q[1];
rz(-1.5953976) q[1];
sx q[1];
rz(0.66397159) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85504882) q[3];
sx q[3];
rz(-2.2576393) q[3];
sx q[3];
rz(-3.1058725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17994443) q[2];
sx q[2];
rz(-2.0521995) q[2];
sx q[2];
rz(1.3570448) q[2];
rz(2.0876136) q[3];
sx q[3];
rz(-1.5318003) q[3];
sx q[3];
rz(2.0747144) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8992952) q[0];
sx q[0];
rz(-0.59277788) q[0];
sx q[0];
rz(1.5574667) q[0];
rz(-1.2552931) q[1];
sx q[1];
rz(-0.86304945) q[1];
sx q[1];
rz(0.042162808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164207) q[0];
sx q[0];
rz(-0.76043425) q[0];
sx q[0];
rz(0.50559931) q[0];
rz(-pi) q[1];
rz(-0.73017786) q[2];
sx q[2];
rz(-0.71766254) q[2];
sx q[2];
rz(-2.4210986) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0706341) q[1];
sx q[1];
rz(-1.8629304) q[1];
sx q[1];
rz(-0.16865428) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8852127) q[3];
sx q[3];
rz(-2.5605128) q[3];
sx q[3];
rz(0.76900859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2178847) q[2];
sx q[2];
rz(-1.8929241) q[2];
sx q[2];
rz(-2.5200747) q[2];
rz(0.38226852) q[3];
sx q[3];
rz(-1.9999802) q[3];
sx q[3];
rz(-2.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650448) q[0];
sx q[0];
rz(-1.5060197) q[0];
sx q[0];
rz(-1.7682834) q[0];
rz(-0.37257591) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(3.0212044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8316702) q[0];
sx q[0];
rz(-0.86747456) q[0];
sx q[0];
rz(2.4256718) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8502813) q[2];
sx q[2];
rz(-1.5320383) q[2];
sx q[2];
rz(-1.9732158) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87625831) q[1];
sx q[1];
rz(-1.8855394) q[1];
sx q[1];
rz(-0.77074428) q[1];
x q[2];
rz(2.1028768) q[3];
sx q[3];
rz(-1.8076583) q[3];
sx q[3];
rz(0.069687931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4176214) q[2];
sx q[2];
rz(-1.3105023) q[2];
sx q[2];
rz(1.854678) q[2];
rz(-0.78479615) q[3];
sx q[3];
rz(-1.1104106) q[3];
sx q[3];
rz(-0.77814656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50412905) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(0.24566393) q[0];
rz(0.0044814666) q[1];
sx q[1];
rz(-2.5878398) q[1];
sx q[1];
rz(3.0779238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5144062) q[0];
sx q[0];
rz(-2.2133491) q[0];
sx q[0];
rz(1.8802934) q[0];
rz(2.0374299) q[2];
sx q[2];
rz(-1.2793102) q[2];
sx q[2];
rz(-1.9574036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88118756) q[1];
sx q[1];
rz(-2.4266325) q[1];
sx q[1];
rz(-1.2843389) q[1];
x q[2];
rz(2.7400209) q[3];
sx q[3];
rz(-2.662559) q[3];
sx q[3];
rz(3.1216321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29023805) q[2];
sx q[2];
rz(-1.3685702) q[2];
sx q[2];
rz(1.9264539) q[2];
rz(0.75567192) q[3];
sx q[3];
rz(-2.1561421) q[3];
sx q[3];
rz(0.22172609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3168685) q[0];
sx q[0];
rz(-1.8727973) q[0];
sx q[0];
rz(0.49279898) q[0];
rz(0.97152501) q[1];
sx q[1];
rz(-1.8172455) q[1];
sx q[1];
rz(-2.9528101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072666) q[0];
sx q[0];
rz(-1.8456689) q[0];
sx q[0];
rz(1.7455533) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9813779) q[2];
sx q[2];
rz(-3.1268178) q[2];
sx q[2];
rz(-2.229634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35313126) q[1];
sx q[1];
rz(-1.6470688) q[1];
sx q[1];
rz(1.0623054) q[1];
rz(-2.4972992) q[3];
sx q[3];
rz(-1.3683934) q[3];
sx q[3];
rz(2.6930893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2762974) q[2];
sx q[2];
rz(-1.5027639) q[2];
sx q[2];
rz(3.0740645) q[2];
rz(-1.6826132) q[3];
sx q[3];
rz(-2.429481) q[3];
sx q[3];
rz(2.0290831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35906288) q[0];
sx q[0];
rz(-1.1787865) q[0];
sx q[0];
rz(-1.68574) q[0];
rz(2.9071232) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(-0.04105982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11066786) q[0];
sx q[0];
rz(-2.2718856) q[0];
sx q[0];
rz(-1.5337471) q[0];
rz(-pi) q[1];
rz(1.7975989) q[2];
sx q[2];
rz(-2.5580346) q[2];
sx q[2];
rz(2.9160519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.252255) q[1];
sx q[1];
rz(-2.4790211) q[1];
sx q[1];
rz(-2.322648) q[1];
rz(-pi) q[2];
rz(1.4047296) q[3];
sx q[3];
rz(-2.4061446) q[3];
sx q[3];
rz(-2.0388132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91512758) q[2];
sx q[2];
rz(-2.3751986) q[2];
sx q[2];
rz(2.8506632) q[2];
rz(-0.21864024) q[3];
sx q[3];
rz(-2.4509957) q[3];
sx q[3];
rz(-0.87029988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17756322) q[0];
sx q[0];
rz(-2.18631) q[0];
sx q[0];
rz(-2.5780504) q[0];
rz(-0.82849416) q[1];
sx q[1];
rz(-2.4528613) q[1];
sx q[1];
rz(-2.4019737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1391382) q[0];
sx q[0];
rz(-1.0268624) q[0];
sx q[0];
rz(-0.55940658) q[0];
x q[1];
rz(0.17156667) q[2];
sx q[2];
rz(-1.0597777) q[2];
sx q[2];
rz(1.9009862) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7259638) q[1];
sx q[1];
rz(-2.0705372) q[1];
sx q[1];
rz(2.8377297) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4007934) q[3];
sx q[3];
rz(-2.1925266) q[3];
sx q[3];
rz(0.19843693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9699817) q[2];
sx q[2];
rz(-2.3889399) q[2];
sx q[2];
rz(-2.8438711) q[2];
rz(0.026084829) q[3];
sx q[3];
rz(-1.9484768) q[3];
sx q[3];
rz(0.74371663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17224269) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(1.7283537) q[0];
rz(-2.6353432) q[1];
sx q[1];
rz(-1.0289611) q[1];
sx q[1];
rz(1.9821573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22102872) q[0];
sx q[0];
rz(-2.8205296) q[0];
sx q[0];
rz(0.32717021) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0721092) q[2];
sx q[2];
rz(-0.58260703) q[2];
sx q[2];
rz(0.16361388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7680507) q[1];
sx q[1];
rz(-0.53540472) q[1];
sx q[1];
rz(0.47799503) q[1];
rz(-pi) q[2];
rz(1.4425464) q[3];
sx q[3];
rz(-0.82144605) q[3];
sx q[3];
rz(-1.8260223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1737698) q[2];
sx q[2];
rz(-0.92542595) q[2];
sx q[2];
rz(-1.3646763) q[2];
rz(2.3260498) q[3];
sx q[3];
rz(-1.848685) q[3];
sx q[3];
rz(1.3968141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5705465) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(1.1676769) q[0];
rz(2.9171464) q[1];
sx q[1];
rz(-2.3851604) q[1];
sx q[1];
rz(0.44463739) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8380008) q[0];
sx q[0];
rz(-1.5193257) q[0];
sx q[0];
rz(-3.0862242) q[0];
rz(-pi) q[1];
rz(-2.3063956) q[2];
sx q[2];
rz(-0.85934256) q[2];
sx q[2];
rz(-0.85853133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9079355) q[1];
sx q[1];
rz(-0.3668712) q[1];
sx q[1];
rz(1.1790102) q[1];
x q[2];
rz(-1.7135213) q[3];
sx q[3];
rz(-1.1253353) q[3];
sx q[3];
rz(0.023825432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4148407) q[2];
sx q[2];
rz(-1.8291992) q[2];
sx q[2];
rz(-1.7337588) q[2];
rz(1.6057181) q[3];
sx q[3];
rz(-1.9460257) q[3];
sx q[3];
rz(-0.87299743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51215148) q[0];
sx q[0];
rz(-2.3245071) q[0];
sx q[0];
rz(-1.0907115) q[0];
rz(2.0953983) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(-2.2549021) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8532904) q[0];
sx q[0];
rz(-1.9294943) q[0];
sx q[0];
rz(-2.7526998) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9891631) q[2];
sx q[2];
rz(-2.5083087) q[2];
sx q[2];
rz(-0.12073244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1266409) q[1];
sx q[1];
rz(-2.0021582) q[1];
sx q[1];
rz(1.5853314) q[1];
x q[2];
rz(-2.7575023) q[3];
sx q[3];
rz(-0.94058296) q[3];
sx q[3];
rz(1.5767136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0118759) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(1.8941619) q[2];
rz(-0.85793197) q[3];
sx q[3];
rz(-1.6777638) q[3];
sx q[3];
rz(-1.5425382) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3457376) q[0];
sx q[0];
rz(-1.7248187) q[0];
sx q[0];
rz(3.0378573) q[0];
rz(0.70032447) q[1];
sx q[1];
rz(-1.843597) q[1];
sx q[1];
rz(1.2784943) q[1];
rz(-2.5896435) q[2];
sx q[2];
rz(-0.23005869) q[2];
sx q[2];
rz(2.1800359) q[2];
rz(-0.92209254) q[3];
sx q[3];
rz(-2.5917883) q[3];
sx q[3];
rz(-2.5592309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
