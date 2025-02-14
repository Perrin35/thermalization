OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(-2.804629) q[0];
rz(0.26951867) q[1];
sx q[1];
rz(-2.1989006) q[1];
sx q[1];
rz(1.2595133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058858697) q[0];
sx q[0];
rz(-1.0089416) q[0];
sx q[0];
rz(0.64942067) q[0];
rz(-0.23064166) q[2];
sx q[2];
rz(-2.3911124) q[2];
sx q[2];
rz(-2.7265158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9971157) q[1];
sx q[1];
rz(-1.0273178) q[1];
sx q[1];
rz(-1.687084) q[1];
rz(2.1416177) q[3];
sx q[3];
rz(-0.68203841) q[3];
sx q[3];
rz(1.2004971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17244615) q[2];
sx q[2];
rz(-1.5772051) q[2];
sx q[2];
rz(-1.2338314) q[2];
rz(1.5305758) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(-0.56510258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3259786) q[0];
sx q[0];
rz(-0.9650721) q[0];
sx q[0];
rz(-2.8976029) q[0];
rz(-1.1583534) q[1];
sx q[1];
rz(-2.127425) q[1];
sx q[1];
rz(-0.44833952) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2650484) q[0];
sx q[0];
rz(-2.6175313) q[0];
sx q[0];
rz(0.35851176) q[0];
rz(0.93083207) q[2];
sx q[2];
rz(-0.85083598) q[2];
sx q[2];
rz(-2.623327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9354466) q[1];
sx q[1];
rz(-2.0191158) q[1];
sx q[1];
rz(2.1765393) q[1];
rz(-pi) q[2];
rz(-0.48359343) q[3];
sx q[3];
rz(-2.0663056) q[3];
sx q[3];
rz(2.7413975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9490732) q[2];
sx q[2];
rz(-3.101888) q[2];
sx q[2];
rz(2.330759) q[2];
rz(-0.0090946322) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(0.31753376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119165) q[0];
sx q[0];
rz(-1.3443953) q[0];
sx q[0];
rz(-2.1938531) q[0];
rz(0.89836994) q[1];
sx q[1];
rz(-0.71304524) q[1];
sx q[1];
rz(-2.9352442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2284828) q[0];
sx q[0];
rz(-2.3867199) q[0];
sx q[0];
rz(-1.55964) q[0];
rz(-pi) q[1];
rz(2.3487665) q[2];
sx q[2];
rz(-1.8603051) q[2];
sx q[2];
rz(0.0057366554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5134597) q[1];
sx q[1];
rz(-1.1185137) q[1];
sx q[1];
rz(-2.9123405) q[1];
rz(-0.32636362) q[3];
sx q[3];
rz(-1.5263399) q[3];
sx q[3];
rz(1.3134559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72948939) q[2];
sx q[2];
rz(-2.1495843) q[2];
sx q[2];
rz(-0.98423973) q[2];
rz(-2.6442243) q[3];
sx q[3];
rz(-1.8822742) q[3];
sx q[3];
rz(-0.74500144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618149) q[0];
sx q[0];
rz(-0.091876939) q[0];
sx q[0];
rz(-2.9352557) q[0];
rz(-1.4641209) q[1];
sx q[1];
rz(-2.3787777) q[1];
sx q[1];
rz(-0.62613553) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3040309) q[0];
sx q[0];
rz(-1.921424) q[0];
sx q[0];
rz(2.3239345) q[0];
rz(-pi) q[1];
x q[1];
rz(3.038806) q[2];
sx q[2];
rz(-1.572552) q[2];
sx q[2];
rz(0.87550113) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5650249) q[1];
sx q[1];
rz(-2.3901689) q[1];
sx q[1];
rz(-2.4016301) q[1];
rz(-pi) q[2];
x q[2];
rz(3.042661) q[3];
sx q[3];
rz(-1.1857934) q[3];
sx q[3];
rz(1.0855293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3760486) q[2];
sx q[2];
rz(-2.3638793) q[2];
sx q[2];
rz(0.98975873) q[2];
rz(-2.1073714) q[3];
sx q[3];
rz(-1.9725622) q[3];
sx q[3];
rz(0.52556747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4393716) q[0];
sx q[0];
rz(-1.419743) q[0];
sx q[0];
rz(-1.1628344) q[0];
rz(-2.974466) q[1];
sx q[1];
rz(-1.5252599) q[1];
sx q[1];
rz(-2.4662245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57608312) q[0];
sx q[0];
rz(-1.352549) q[0];
sx q[0];
rz(-0.670114) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.295601) q[2];
sx q[2];
rz(-0.48946871) q[2];
sx q[2];
rz(-0.54014017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4247508) q[1];
sx q[1];
rz(-1.9159686) q[1];
sx q[1];
rz(1.1139681) q[1];
x q[2];
rz(2.1618103) q[3];
sx q[3];
rz(-1.3434935) q[3];
sx q[3];
rz(2.1066372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2165788) q[2];
sx q[2];
rz(-1.4224956) q[2];
sx q[2];
rz(2.6340346) q[2];
rz(-0.01928586) q[3];
sx q[3];
rz(-2.3465893) q[3];
sx q[3];
rz(1.9222586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5386706) q[0];
sx q[0];
rz(-0.62748533) q[0];
sx q[0];
rz(-0.70166171) q[0];
rz(0.64274669) q[1];
sx q[1];
rz(-1.9822491) q[1];
sx q[1];
rz(1.3210375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82133991) q[0];
sx q[0];
rz(-2.5139315) q[0];
sx q[0];
rz(-1.5492803) q[0];
x q[1];
rz(-2.7847544) q[2];
sx q[2];
rz(-0.90274397) q[2];
sx q[2];
rz(-0.84260637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63032645) q[1];
sx q[1];
rz(-2.0980586) q[1];
sx q[1];
rz(0.43095891) q[1];
rz(-2.3307468) q[3];
sx q[3];
rz(-1.7827991) q[3];
sx q[3];
rz(-0.9190587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.295149) q[2];
sx q[2];
rz(-0.48131338) q[2];
sx q[2];
rz(-2.5018137) q[2];
rz(-2.4844737) q[3];
sx q[3];
rz(-1.6491456) q[3];
sx q[3];
rz(-2.8583728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.101864) q[0];
sx q[0];
rz(-1.855408) q[0];
sx q[0];
rz(-0.78045994) q[0];
rz(-1.9038433) q[1];
sx q[1];
rz(-1.8433808) q[1];
sx q[1];
rz(-0.16403988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17724161) q[0];
sx q[0];
rz(-0.57946396) q[0];
sx q[0];
rz(2.01418) q[0];
rz(-2.7299147) q[2];
sx q[2];
rz(-1.9814166) q[2];
sx q[2];
rz(-3.0268998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32719993) q[1];
sx q[1];
rz(-2.3770824) q[1];
sx q[1];
rz(2.3491377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79188723) q[3];
sx q[3];
rz(-1.2372176) q[3];
sx q[3];
rz(1.4871979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8728472) q[2];
sx q[2];
rz(-1.9445323) q[2];
sx q[2];
rz(-2.2609113) q[2];
rz(1.4745447) q[3];
sx q[3];
rz(-1.0677974) q[3];
sx q[3];
rz(1.4145981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6678161) q[0];
sx q[0];
rz(-2.9528463) q[0];
sx q[0];
rz(1.1291946) q[0];
rz(0.95343626) q[1];
sx q[1];
rz(-0.9402746) q[1];
sx q[1];
rz(-0.58852351) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584483) q[0];
sx q[0];
rz(-2.5059688) q[0];
sx q[0];
rz(1.0516214) q[0];
x q[1];
rz(1.2382201) q[2];
sx q[2];
rz(-1.7875449) q[2];
sx q[2];
rz(-2.7970566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.037884574) q[1];
sx q[1];
rz(-2.3988535) q[1];
sx q[1];
rz(2.6602416) q[1];
rz(-pi) q[2];
rz(1.3422119) q[3];
sx q[3];
rz(-0.22998631) q[3];
sx q[3];
rz(-3.1222469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15527655) q[2];
sx q[2];
rz(-2.1337815) q[2];
sx q[2];
rz(1.0922095) q[2];
rz(2.3327667) q[3];
sx q[3];
rz(-2.0598965) q[3];
sx q[3];
rz(1.697418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5876193) q[0];
sx q[0];
rz(-1.076979) q[0];
sx q[0];
rz(0.67163604) q[0];
rz(0.49297586) q[1];
sx q[1];
rz(-2.2608345) q[1];
sx q[1];
rz(-1.2395476) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9818891) q[0];
sx q[0];
rz(-2.9043439) q[0];
sx q[0];
rz(-1.5818198) q[0];
x q[1];
rz(-2.9443594) q[2];
sx q[2];
rz(-1.5610362) q[2];
sx q[2];
rz(2.724144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0555435) q[1];
sx q[1];
rz(-0.99850875) q[1];
sx q[1];
rz(-1.0466762) q[1];
x q[2];
rz(-2.3726241) q[3];
sx q[3];
rz(-1.4533224) q[3];
sx q[3];
rz(-2.3905135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12941831) q[2];
sx q[2];
rz(-1.6678026) q[2];
sx q[2];
rz(2.5950281) q[2];
rz(0.92285815) q[3];
sx q[3];
rz(-0.62896362) q[3];
sx q[3];
rz(-1.9950689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65761143) q[0];
sx q[0];
rz(-0.46009362) q[0];
sx q[0];
rz(2.7125603) q[0];
rz(1.8589164) q[1];
sx q[1];
rz(-1.7762643) q[1];
sx q[1];
rz(-0.081258953) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51516303) q[0];
sx q[0];
rz(-1.7395955) q[0];
sx q[0];
rz(-1.4909046) q[0];
rz(-pi) q[1];
rz(0.0325412) q[2];
sx q[2];
rz(-2.0306808) q[2];
sx q[2];
rz(-2.3655917) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1926905) q[1];
sx q[1];
rz(-1.16445) q[1];
sx q[1];
rz(-2.3188616) q[1];
rz(-pi) q[2];
rz(-2.6762415) q[3];
sx q[3];
rz(-1.1123667) q[3];
sx q[3];
rz(-2.0260594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73605865) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(-2.0299358) q[2];
rz(-1.9723655) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(0.14298239) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0040759857) q[0];
sx q[0];
rz(-0.62340323) q[0];
sx q[0];
rz(2.2829983) q[0];
rz(2.2463592) q[1];
sx q[1];
rz(-1.8276855) q[1];
sx q[1];
rz(-3.102416) q[1];
rz(-1.8276385) q[2];
sx q[2];
rz(-1.5821291) q[2];
sx q[2];
rz(-0.24072992) q[2];
rz(-1.7029892) q[3];
sx q[3];
rz(-0.93441603) q[3];
sx q[3];
rz(-0.84032755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
