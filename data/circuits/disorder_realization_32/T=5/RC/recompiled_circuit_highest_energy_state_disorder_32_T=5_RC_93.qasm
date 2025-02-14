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
rz(-2.3218396) q[1];
sx q[1];
rz(4.0263962) q[1];
sx q[1];
rz(10.560992) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3216517) q[0];
sx q[0];
rz(-2.5188832) q[0];
sx q[0];
rz(1.4371928) q[0];
rz(-2.4769463) q[2];
sx q[2];
rz(-2.342939) q[2];
sx q[2];
rz(-2.2448774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6032049) q[1];
sx q[1];
rz(-0.520272) q[1];
sx q[1];
rz(0.59641312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3435279) q[3];
sx q[3];
rz(-1.1385001) q[3];
sx q[3];
rz(-2.5903518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39584407) q[2];
sx q[2];
rz(-2.6555847) q[2];
sx q[2];
rz(2.8493472) q[2];
rz(-2.5203943) q[3];
sx q[3];
rz(-1.0750333) q[3];
sx q[3];
rz(-0.21875374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.56276965) q[0];
sx q[0];
rz(-2.0243473) q[0];
sx q[0];
rz(-3.0551832) q[0];
rz(0.39331618) q[1];
sx q[1];
rz(-1.4540648) q[1];
sx q[1];
rz(-1.6578065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1817966) q[0];
sx q[0];
rz(-2.2364669) q[0];
sx q[0];
rz(-0.29924198) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1643402) q[2];
sx q[2];
rz(-2.8033442) q[2];
sx q[2];
rz(-1.8601314) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68951974) q[1];
sx q[1];
rz(-1.7240925) q[1];
sx q[1];
rz(0.91233715) q[1];
rz(-1.3065763) q[3];
sx q[3];
rz(-2.1253028) q[3];
sx q[3];
rz(-2.4305024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1958486) q[2];
sx q[2];
rz(-1.5429147) q[2];
sx q[2];
rz(-3.1405385) q[2];
rz(-0.70683181) q[3];
sx q[3];
rz(-0.89444923) q[3];
sx q[3];
rz(2.8640174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23866776) q[0];
sx q[0];
rz(-2.213573) q[0];
sx q[0];
rz(2.6631885) q[0];
rz(-0.84486419) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(2.1873651) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3197495) q[0];
sx q[0];
rz(-1.1992393) q[0];
sx q[0];
rz(-0.48858924) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61561959) q[2];
sx q[2];
rz(-2.2055045) q[2];
sx q[2];
rz(-2.8826432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0597998) q[1];
sx q[1];
rz(-0.89443086) q[1];
sx q[1];
rz(-2.6790819) q[1];
x q[2];
rz(-0.34326633) q[3];
sx q[3];
rz(-2.4606554) q[3];
sx q[3];
rz(-1.9576114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2039631) q[2];
sx q[2];
rz(-0.340168) q[2];
sx q[2];
rz(3.0703239) q[2];
rz(2.1313306) q[3];
sx q[3];
rz(-1.1610169) q[3];
sx q[3];
rz(-2.8205813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45692745) q[0];
sx q[0];
rz(-1.848897) q[0];
sx q[0];
rz(-2.7780823) q[0];
rz(0.89206308) q[1];
sx q[1];
rz(-2.108768) q[1];
sx q[1];
rz(1.8878638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15479408) q[0];
sx q[0];
rz(-2.5240233) q[0];
sx q[0];
rz(-2.2876431) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10102306) q[2];
sx q[2];
rz(-1.2083078) q[2];
sx q[2];
rz(-0.50607888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82361932) q[1];
sx q[1];
rz(-2.8880305) q[1];
sx q[1];
rz(-2.9754354) q[1];
rz(2.4917077) q[3];
sx q[3];
rz(-2.657848) q[3];
sx q[3];
rz(0.080881491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7518647) q[2];
sx q[2];
rz(-2.3675297) q[2];
sx q[2];
rz(-1.6130201) q[2];
rz(0.45016897) q[3];
sx q[3];
rz(-1.3521786) q[3];
sx q[3];
rz(-1.2990797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0171829) q[0];
sx q[0];
rz(-2.9876509) q[0];
sx q[0];
rz(2.2907139) q[0];
rz(-1.7393913) q[1];
sx q[1];
rz(-1.5504928) q[1];
sx q[1];
rz(-2.9891787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35203347) q[0];
sx q[0];
rz(-1.5769582) q[0];
sx q[0];
rz(-1.8101109) q[0];
rz(-pi) q[1];
rz(2.1722069) q[2];
sx q[2];
rz(-1.3297373) q[2];
sx q[2];
rz(-3.0066662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20402699) q[1];
sx q[1];
rz(-1.7364677) q[1];
sx q[1];
rz(-1.8574287) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39800003) q[3];
sx q[3];
rz(-1.9129244) q[3];
sx q[3];
rz(2.0002535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57205814) q[2];
sx q[2];
rz(-2.2172838) q[2];
sx q[2];
rz(2.7657236) q[2];
rz(2.87319) q[3];
sx q[3];
rz(-1.0207876) q[3];
sx q[3];
rz(-2.102803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1035476) q[0];
sx q[0];
rz(-0.91829848) q[0];
sx q[0];
rz(0.12575664) q[0];
rz(-3.0962931) q[1];
sx q[1];
rz(-2.2442975) q[1];
sx q[1];
rz(-0.50517857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99085) q[0];
sx q[0];
rz(-1.8253541) q[0];
sx q[0];
rz(1.6360511) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.061902) q[2];
sx q[2];
rz(-1.4734846) q[2];
sx q[2];
rz(-2.8081913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1152079) q[1];
sx q[1];
rz(-1.4248825) q[1];
sx q[1];
rz(-2.3664631) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.807928) q[3];
sx q[3];
rz(-2.0150321) q[3];
sx q[3];
rz(-2.9315116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3230285) q[2];
sx q[2];
rz(-2.8105141) q[2];
sx q[2];
rz(0.25311145) q[2];
rz(-0.74462914) q[3];
sx q[3];
rz(-1.587845) q[3];
sx q[3];
rz(0.26960069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.700106) q[0];
sx q[0];
rz(-1.8550669) q[0];
sx q[0];
rz(-0.043524608) q[0];
rz(-1.9684017) q[1];
sx q[1];
rz(-1.2930608) q[1];
sx q[1];
rz(-0.80165577) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7907958) q[0];
sx q[0];
rz(-2.9614095) q[0];
sx q[0];
rz(-2.5596746) q[0];
x q[1];
rz(0.72791667) q[2];
sx q[2];
rz(-2.2322828) q[2];
sx q[2];
rz(2.1833411) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3139627) q[1];
sx q[1];
rz(-2.4555742) q[1];
sx q[1];
rz(-1.3242701) q[1];
rz(-1.8160062) q[3];
sx q[3];
rz(-1.86487) q[3];
sx q[3];
rz(0.40286703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39048553) q[2];
sx q[2];
rz(-1.0932873) q[2];
sx q[2];
rz(-2.2170179) q[2];
rz(-3.0604002) q[3];
sx q[3];
rz(-1.743318) q[3];
sx q[3];
rz(-2.413522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28295949) q[0];
sx q[0];
rz(-1.8959683) q[0];
sx q[0];
rz(-2.0661195) q[0];
rz(1.3772427) q[1];
sx q[1];
rz(-1.7694764) q[1];
sx q[1];
rz(-2.4864054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7542962) q[0];
sx q[0];
rz(-1.5167902) q[0];
sx q[0];
rz(-1.61458) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29187188) q[2];
sx q[2];
rz(-0.8093738) q[2];
sx q[2];
rz(0.82631095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1121443) q[1];
sx q[1];
rz(-0.69975425) q[1];
sx q[1];
rz(-3.0359731) q[1];
rz(1.4873679) q[3];
sx q[3];
rz(-2.8018824) q[3];
sx q[3];
rz(1.8172906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7738771) q[2];
sx q[2];
rz(-2.750062) q[2];
sx q[2];
rz(1.1565304) q[2];
rz(-1.997442) q[3];
sx q[3];
rz(-0.63026989) q[3];
sx q[3];
rz(1.8088079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82041204) q[0];
sx q[0];
rz(-0.92083609) q[0];
sx q[0];
rz(-1.9658827) q[0];
rz(1.3395478) q[1];
sx q[1];
rz(-2.1422155) q[1];
sx q[1];
rz(0.40750009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20243199) q[0];
sx q[0];
rz(-0.5078041) q[0];
sx q[0];
rz(2.5277471) q[0];
x q[1];
rz(2.6596498) q[2];
sx q[2];
rz(-0.53496581) q[2];
sx q[2];
rz(-0.15967655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6432156) q[1];
sx q[1];
rz(-0.98340323) q[1];
sx q[1];
rz(0.024557928) q[1];
rz(-pi) q[2];
rz(-0.27176933) q[3];
sx q[3];
rz(-1.5021694) q[3];
sx q[3];
rz(2.8420699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0672062) q[2];
sx q[2];
rz(-1.1585453) q[2];
sx q[2];
rz(-1.0254394) q[2];
rz(0.49753672) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(2.735125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24791726) q[0];
sx q[0];
rz(-1.9837288) q[0];
sx q[0];
rz(-1.3854618) q[0];
rz(-2.0939743) q[1];
sx q[1];
rz(-1.3448998) q[1];
sx q[1];
rz(-2.1687543) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6167283) q[0];
sx q[0];
rz(-1.4848733) q[0];
sx q[0];
rz(2.1154386) q[0];
x q[1];
rz(-0.49373105) q[2];
sx q[2];
rz(-2.0539453) q[2];
sx q[2];
rz(-1.2497304) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23215881) q[1];
sx q[1];
rz(-1.9240018) q[1];
sx q[1];
rz(-1.722419) q[1];
x q[2];
rz(-0.38145251) q[3];
sx q[3];
rz(-0.90899639) q[3];
sx q[3];
rz(-1.0009777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4086548) q[2];
sx q[2];
rz(-2.813377) q[2];
sx q[2];
rz(1.1019863) q[2];
rz(-1.7105626) q[3];
sx q[3];
rz(-2.3165063) q[3];
sx q[3];
rz(1.2380606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86193209) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(-2.4109789) q[1];
sx q[1];
rz(-2.5082671) q[1];
sx q[1];
rz(-2.3309753) q[1];
rz(-2.7027086) q[2];
sx q[2];
rz(-0.65541291) q[2];
sx q[2];
rz(2.4272774) q[2];
rz(2.4875658) q[3];
sx q[3];
rz(-0.43523052) q[3];
sx q[3];
rz(-0.38149618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
