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
rz(3.1130528) q[0];
sx q[0];
rz(-0.8172577) q[0];
sx q[0];
rz(0.28491268) q[0];
rz(0.81975308) q[1];
sx q[1];
rz(-0.88480359) q[1];
sx q[1];
rz(2.0053782) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2819913) q[0];
sx q[0];
rz(-1.4930269) q[0];
sx q[0];
rz(-2.1892709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66464632) q[2];
sx q[2];
rz(-0.79865361) q[2];
sx q[2];
rz(0.89671521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49992554) q[1];
sx q[1];
rz(-1.8537775) q[1];
sx q[1];
rz(2.6989515) q[1];
rz(-pi) q[2];
rz(1.3435279) q[3];
sx q[3];
rz(-1.1385001) q[3];
sx q[3];
rz(-0.55124084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39584407) q[2];
sx q[2];
rz(-0.48600799) q[2];
sx q[2];
rz(-2.8493472) q[2];
rz(2.5203943) q[3];
sx q[3];
rz(-1.0750333) q[3];
sx q[3];
rz(0.21875374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.578823) q[0];
sx q[0];
rz(-1.1172453) q[0];
sx q[0];
rz(-3.0551832) q[0];
rz(-0.39331618) q[1];
sx q[1];
rz(-1.4540648) q[1];
sx q[1];
rz(1.6578065) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.718865) q[0];
sx q[0];
rz(-1.3368092) q[0];
sx q[0];
rz(-0.88293332) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8833161) q[2];
sx q[2];
rz(-1.7023689) q[2];
sx q[2];
rz(2.4665592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4520729) q[1];
sx q[1];
rz(-1.7240925) q[1];
sx q[1];
rz(2.2292555) q[1];
rz(-0.57043907) q[3];
sx q[3];
rz(-1.794687) q[3];
sx q[3];
rz(-0.71820948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.94574404) q[2];
sx q[2];
rz(-1.5986779) q[2];
sx q[2];
rz(3.1405385) q[2];
rz(-0.70683181) q[3];
sx q[3];
rz(-0.89444923) q[3];
sx q[3];
rz(2.8640174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23866776) q[0];
sx q[0];
rz(-2.213573) q[0];
sx q[0];
rz(-0.4784041) q[0];
rz(0.84486419) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(0.9542276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8218431) q[0];
sx q[0];
rz(-1.9423534) q[0];
sx q[0];
rz(2.6530034) q[0];
rz(-pi) q[1];
rz(-0.61561959) q[2];
sx q[2];
rz(-2.2055045) q[2];
sx q[2];
rz(0.25894946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79151691) q[1];
sx q[1];
rz(-1.215394) q[1];
sx q[1];
rz(-2.3018964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65173463) q[3];
sx q[3];
rz(-1.7842891) q[3];
sx q[3];
rz(-0.65769559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93762952) q[2];
sx q[2];
rz(-0.340168) q[2];
sx q[2];
rz(-0.071268737) q[2];
rz(2.1313306) q[3];
sx q[3];
rz(-1.1610169) q[3];
sx q[3];
rz(0.32101139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6846652) q[0];
sx q[0];
rz(-1.848897) q[0];
sx q[0];
rz(0.36351031) q[0];
rz(-0.89206308) q[1];
sx q[1];
rz(-2.108768) q[1];
sx q[1];
rz(-1.8878638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79824582) q[0];
sx q[0];
rz(-1.1805184) q[0];
sx q[0];
rz(1.0791995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3108861) q[2];
sx q[2];
rz(-0.37570243) q[2];
sx q[2];
rz(2.9139522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82361932) q[1];
sx q[1];
rz(-0.25356217) q[1];
sx q[1];
rz(-2.9754354) q[1];
rz(0.39616743) q[3];
sx q[3];
rz(-1.2855144) q[3];
sx q[3];
rz(1.0594291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38972798) q[2];
sx q[2];
rz(-0.77406293) q[2];
sx q[2];
rz(1.6130201) q[2];
rz(-2.6914237) q[3];
sx q[3];
rz(-1.7894141) q[3];
sx q[3];
rz(-1.842513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12440974) q[0];
sx q[0];
rz(-2.9876509) q[0];
sx q[0];
rz(-2.2907139) q[0];
rz(-1.7393913) q[1];
sx q[1];
rz(-1.5910999) q[1];
sx q[1];
rz(2.9891787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440118) q[0];
sx q[0];
rz(-0.23939238) q[0];
sx q[0];
rz(-1.5448065) q[0];
rz(-pi) q[1];
rz(-1.1609116) q[2];
sx q[2];
rz(-0.64233795) q[2];
sx q[2];
rz(-1.1010686) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7262532) q[1];
sx q[1];
rz(-1.2881973) q[1];
sx q[1];
rz(0.17258172) q[1];
x q[2];
rz(1.9394623) q[3];
sx q[3];
rz(-1.1970425) q[3];
sx q[3];
rz(2.8522648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5695345) q[2];
sx q[2];
rz(-2.2172838) q[2];
sx q[2];
rz(0.37586907) q[2];
rz(-2.87319) q[3];
sx q[3];
rz(-1.0207876) q[3];
sx q[3];
rz(2.102803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.038045) q[0];
sx q[0];
rz(-0.91829848) q[0];
sx q[0];
rz(0.12575664) q[0];
rz(-3.0962931) q[1];
sx q[1];
rz(-2.2442975) q[1];
sx q[1];
rz(2.6364141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4046458) q[0];
sx q[0];
rz(-2.8789799) q[0];
sx q[0];
rz(-0.24554952) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11025496) q[2];
sx q[2];
rz(-1.0822191) q[2];
sx q[2];
rz(-1.9561121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40292254) q[1];
sx q[1];
rz(-0.80602491) q[1];
sx q[1];
rz(-1.3678985) q[1];
rz(-pi) q[2];
rz(2.1734851) q[3];
sx q[3];
rz(-0.54881964) q[3];
sx q[3];
rz(0.88879648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81856412) q[2];
sx q[2];
rz(-2.8105141) q[2];
sx q[2];
rz(0.25311145) q[2];
rz(2.3969635) q[3];
sx q[3];
rz(-1.5537477) q[3];
sx q[3];
rz(2.871992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.44148663) q[0];
sx q[0];
rz(-1.2865257) q[0];
sx q[0];
rz(3.098068) q[0];
rz(-1.9684017) q[1];
sx q[1];
rz(-1.8485319) q[1];
sx q[1];
rz(0.80165577) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2013595) q[0];
sx q[0];
rz(-1.721075) q[0];
sx q[0];
rz(1.6705832) q[0];
x q[1];
rz(-0.76446587) q[2];
sx q[2];
rz(-1.0180961) q[2];
sx q[2];
rz(-2.0281731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55085631) q[1];
sx q[1];
rz(-1.4155861) q[1];
sx q[1];
rz(2.2418145) q[1];
rz(-pi) q[2];
rz(-0.67569895) q[3];
sx q[3];
rz(-2.7609918) q[3];
sx q[3];
rz(0.30932793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7511071) q[2];
sx q[2];
rz(-1.0932873) q[2];
sx q[2];
rz(2.2170179) q[2];
rz(3.0604002) q[3];
sx q[3];
rz(-1.3982747) q[3];
sx q[3];
rz(0.72807062) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28295949) q[0];
sx q[0];
rz(-1.8959683) q[0];
sx q[0];
rz(2.0661195) q[0];
rz(-1.7643499) q[1];
sx q[1];
rz(-1.7694764) q[1];
sx q[1];
rz(-2.4864054) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1811349) q[0];
sx q[0];
rz(-1.5270766) q[0];
sx q[0];
rz(3.0875348) q[0];
rz(-pi) q[1];
rz(0.78777625) q[2];
sx q[2];
rz(-1.3609741) q[2];
sx q[2];
rz(2.1927046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1671518) q[1];
sx q[1];
rz(-0.87572423) q[1];
sx q[1];
rz(1.6593169) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4873679) q[3];
sx q[3];
rz(-2.8018824) q[3];
sx q[3];
rz(1.8172906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7738771) q[2];
sx q[2];
rz(-2.750062) q[2];
sx q[2];
rz(1.9850622) q[2];
rz(-1.997442) q[3];
sx q[3];
rz(-2.5113228) q[3];
sx q[3];
rz(1.3327848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82041204) q[0];
sx q[0];
rz(-0.92083609) q[0];
sx q[0];
rz(1.9658827) q[0];
rz(1.8020449) q[1];
sx q[1];
rz(-2.1422155) q[1];
sx q[1];
rz(2.7340926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88107318) q[0];
sx q[0];
rz(-1.1620191) q[0];
sx q[0];
rz(-1.2606032) q[0];
rz(-pi) q[1];
rz(1.3027329) q[2];
sx q[2];
rz(-1.1020793) q[2];
sx q[2];
rz(2.7549636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.54267) q[1];
sx q[1];
rz(-2.5537468) q[1];
sx q[1];
rz(-1.6076615) q[1];
rz(-pi) q[2];
rz(-0.27176933) q[3];
sx q[3];
rz(-1.5021694) q[3];
sx q[3];
rz(2.8420699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0672062) q[2];
sx q[2];
rz(-1.1585453) q[2];
sx q[2];
rz(1.0254394) q[2];
rz(-2.6440559) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(2.735125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936754) q[0];
sx q[0];
rz(-1.9837288) q[0];
sx q[0];
rz(-1.7561308) q[0];
rz(-2.0939743) q[1];
sx q[1];
rz(-1.3448998) q[1];
sx q[1];
rz(0.97283831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0437229) q[0];
sx q[0];
rz(-2.1132054) q[0];
sx q[0];
rz(-0.10036758) q[0];
rz(0.83613427) q[2];
sx q[2];
rz(-0.67648602) q[2];
sx q[2];
rz(-2.1084146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23215881) q[1];
sx q[1];
rz(-1.9240018) q[1];
sx q[1];
rz(1.722419) q[1];
rz(-pi) q[2];
rz(-1.1249969) q[3];
sx q[3];
rz(-0.74927038) q[3];
sx q[3];
rz(-0.42271915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73293781) q[2];
sx q[2];
rz(-2.813377) q[2];
sx q[2];
rz(2.0396063) q[2];
rz(-1.7105626) q[3];
sx q[3];
rz(-0.82508636) q[3];
sx q[3];
rz(1.903532) q[3];
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
rz(pi/2) q[0];
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
rz(2.2796606) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(0.73061371) q[1];
sx q[1];
rz(-2.5082671) q[1];
sx q[1];
rz(-2.3309753) q[1];
rz(-1.8865449) q[2];
sx q[2];
rz(-0.98636711) q[2];
sx q[2];
rz(-1.2489088) q[2];
rz(2.788078) q[3];
sx q[3];
rz(-1.3113889) q[3];
sx q[3];
rz(1.7967381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
