OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(-0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44088988) q[0];
sx q[0];
rz(-2.7452677) q[0];
sx q[0];
rz(1.2896145) q[0];
x q[1];
rz(2.2719703) q[2];
sx q[2];
rz(-0.50719559) q[2];
sx q[2];
rz(-1.6091572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5296386) q[1];
sx q[1];
rz(-2.0679592) q[1];
sx q[1];
rz(1.0477209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3832541) q[3];
sx q[3];
rz(-1.7539277) q[3];
sx q[3];
rz(1.6299712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-0.83591998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0854557) q[0];
sx q[0];
rz(-1.5526999) q[0];
sx q[0];
rz(1.4319112) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95401986) q[2];
sx q[2];
rz(-2.4980133) q[2];
sx q[2];
rz(1.2611024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0554725) q[1];
sx q[1];
rz(-1.8430084) q[1];
sx q[1];
rz(-0.91113669) q[1];
rz(-pi) q[2];
rz(2.4017176) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(-1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-1.0773766) q[0];
sx q[0];
rz(2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3041829) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(-3.0512179) q[0];
rz(1.978546) q[2];
sx q[2];
rz(-1.7951269) q[2];
sx q[2];
rz(3.0293904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2512974) q[1];
sx q[1];
rz(-1.9922868) q[1];
sx q[1];
rz(1.4189659) q[1];
x q[2];
rz(1.6821074) q[3];
sx q[3];
rz(-0.49735945) q[3];
sx q[3];
rz(1.5745844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1221216) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(-2.1602901) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(2.9503126) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0606196) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(-3.0540375) q[0];
rz(0.83696604) q[2];
sx q[2];
rz(-1.3651197) q[2];
sx q[2];
rz(-0.93081805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0335238) q[1];
sx q[1];
rz(-2.3571157) q[1];
sx q[1];
rz(2.9963521) q[1];
x q[2];
rz(0.53381613) q[3];
sx q[3];
rz(-1.8752408) q[3];
sx q[3];
rz(-0.75418762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-2.5849735) q[0];
rz(-3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(0.97250485) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9113377) q[0];
sx q[0];
rz(-1.6382434) q[0];
sx q[0];
rz(-0.10226843) q[0];
rz(0.82053484) q[2];
sx q[2];
rz(-1.3247326) q[2];
sx q[2];
rz(-0.29124242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3481808) q[1];
sx q[1];
rz(-2.9017397) q[1];
sx q[1];
rz(-0.94437771) q[1];
rz(-pi) q[2];
rz(0.036690849) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3187023) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4218629) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.4259889) q[0];
rz(-2.0772207) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(2.7672966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74852809) q[0];
sx q[0];
rz(-0.52951282) q[0];
sx q[0];
rz(0.66324309) q[0];
x q[1];
rz(2.6475545) q[2];
sx q[2];
rz(-1.3400153) q[2];
sx q[2];
rz(-1.5766694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.018505521) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(-1.2612543) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66773325) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(-2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(-1.8315171) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10446564) q[0];
sx q[0];
rz(-1.5195527) q[0];
sx q[0];
rz(-0.38802223) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2398948) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(0.14724018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9359365) q[1];
sx q[1];
rz(-1.2381136) q[1];
sx q[1];
rz(1.3693621) q[1];
x q[2];
rz(-2.5861916) q[3];
sx q[3];
rz(-1.5474833) q[3];
sx q[3];
rz(2.5454552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(-0.84987744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4655315) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(-1.3423052) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2539027) q[2];
sx q[2];
rz(-1.4978834) q[2];
sx q[2];
rz(0.86276744) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1666607) q[1];
sx q[1];
rz(-2.3798675) q[1];
sx q[1];
rz(-1.6208266) q[1];
rz(1.8306584) q[3];
sx q[3];
rz(-1.4514918) q[3];
sx q[3];
rz(-2.3343149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(-2.3855709) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-0.35593629) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-0.40503043) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(-1.8639494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67280806) q[0];
sx q[0];
rz(-2.0622258) q[0];
sx q[0];
rz(1.3093033) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.054504) q[2];
sx q[2];
rz(-2.2030369) q[2];
sx q[2];
rz(1.5916923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1956049) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(3.1006378) q[1];
rz(0.83417474) q[3];
sx q[3];
rz(-1.8636384) q[3];
sx q[3];
rz(-2.1017696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-0.35783106) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(-0.21044883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99104098) q[0];
sx q[0];
rz(-1.7518839) q[0];
sx q[0];
rz(-0.51245706) q[0];
x q[1];
rz(0.16410447) q[2];
sx q[2];
rz(-1.9774984) q[2];
sx q[2];
rz(-1.1526398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3810972) q[1];
sx q[1];
rz(-1.7777182) q[1];
sx q[1];
rz(-0.31462545) q[1];
x q[2];
rz(-0.25037346) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(-1.8680343) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(2.9343228) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(0.81644425) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(-2.0408761) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(2.3002426) q[3];
sx q[3];
rz(-1.8302866) q[3];
sx q[3];
rz(2.4091099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];