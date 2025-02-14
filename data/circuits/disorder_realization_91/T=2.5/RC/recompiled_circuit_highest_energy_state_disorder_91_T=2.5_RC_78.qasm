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
rz(0.36747992) q[0];
sx q[0];
rz(4.9541192) q[0];
sx q[0];
rz(10.911565) q[0];
rz(-1.6410671) q[1];
sx q[1];
rz(-0.60368109) q[1];
sx q[1];
rz(0.85890213) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4444111) q[0];
sx q[0];
rz(-2.1281805) q[0];
sx q[0];
rz(2.0472705) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3043465) q[2];
sx q[2];
rz(-0.7310125) q[2];
sx q[2];
rz(1.3251418) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6851226) q[1];
sx q[1];
rz(-0.89159144) q[1];
sx q[1];
rz(0.14107708) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9263622) q[3];
sx q[3];
rz(-1.269333) q[3];
sx q[3];
rz(1.3380443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19771244) q[2];
sx q[2];
rz(-0.79215017) q[2];
sx q[2];
rz(-2.0749178) q[2];
rz(-0.094430447) q[3];
sx q[3];
rz(-2.5240199) q[3];
sx q[3];
rz(2.7134231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601032) q[0];
sx q[0];
rz(-2.8747989) q[0];
sx q[0];
rz(1.9285404) q[0];
rz(-0.45252291) q[1];
sx q[1];
rz(-2.8892543) q[1];
sx q[1];
rz(0.67950335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247422) q[0];
sx q[0];
rz(-2.8899414) q[0];
sx q[0];
rz(-0.13435039) q[0];
rz(2.8507336) q[2];
sx q[2];
rz(-1.4122987) q[2];
sx q[2];
rz(3.045451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6399709) q[1];
sx q[1];
rz(-2.2414329) q[1];
sx q[1];
rz(2.0790504) q[1];
x q[2];
rz(-0.68778681) q[3];
sx q[3];
rz(-1.4151598) q[3];
sx q[3];
rz(-0.17222541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54000336) q[2];
sx q[2];
rz(-2.3175779) q[2];
sx q[2];
rz(-0.69671112) q[2];
rz(-1.2434897) q[3];
sx q[3];
rz(-1.1949298) q[3];
sx q[3];
rz(0.5881601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.8884647) q[0];
sx q[0];
rz(-1.0364113) q[0];
sx q[0];
rz(0.82365197) q[0];
rz(-0.15692391) q[1];
sx q[1];
rz(-1.7053968) q[1];
sx q[1];
rz(-2.7395111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0500129) q[0];
sx q[0];
rz(-2.7224143) q[0];
sx q[0];
rz(-1.242595) q[0];
x q[1];
rz(1.4060611) q[2];
sx q[2];
rz(-2.852612) q[2];
sx q[2];
rz(1.8584205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7618708) q[1];
sx q[1];
rz(-2.5609697) q[1];
sx q[1];
rz(2.1942744) q[1];
rz(-pi) q[2];
rz(-1.1810494) q[3];
sx q[3];
rz(-1.3661091) q[3];
sx q[3];
rz(-1.8738511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0154401) q[2];
sx q[2];
rz(-2.2045279) q[2];
sx q[2];
rz(0.13119571) q[2];
rz(2.0845856) q[3];
sx q[3];
rz(-1.9481235) q[3];
sx q[3];
rz(-1.2098275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4256725) q[0];
sx q[0];
rz(-1.1765867) q[0];
sx q[0];
rz(-1.7536989) q[0];
rz(-1.362494) q[1];
sx q[1];
rz(-1.250993) q[1];
sx q[1];
rz(-1.2066134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2817474) q[0];
sx q[0];
rz(-2.3424405) q[0];
sx q[0];
rz(-0.65779442) q[0];
rz(0.64995857) q[2];
sx q[2];
rz(-2.0043249) q[2];
sx q[2];
rz(2.1999846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6304508) q[1];
sx q[1];
rz(-1.2690374) q[1];
sx q[1];
rz(1.4901699) q[1];
x q[2];
rz(1.1251757) q[3];
sx q[3];
rz(-0.88332159) q[3];
sx q[3];
rz(-2.4514707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44986192) q[2];
sx q[2];
rz(-1.3351771) q[2];
sx q[2];
rz(-2.195669) q[2];
rz(-0.45390421) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40189633) q[0];
sx q[0];
rz(-1.8296158) q[0];
sx q[0];
rz(0.78347462) q[0];
rz(2.2354194) q[1];
sx q[1];
rz(-0.89007354) q[1];
sx q[1];
rz(0.56437147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2191129) q[0];
sx q[0];
rz(-1.4267491) q[0];
sx q[0];
rz(-0.12024583) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1644092) q[2];
sx q[2];
rz(-1.5157852) q[2];
sx q[2];
rz(-2.0470574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7184938) q[1];
sx q[1];
rz(-2.1717487) q[1];
sx q[1];
rz(-0.81032958) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71471148) q[3];
sx q[3];
rz(-1.6922975) q[3];
sx q[3];
rz(0.35706676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2169317) q[2];
sx q[2];
rz(-1.4035808) q[2];
sx q[2];
rz(3.0972287) q[2];
rz(1.7313322) q[3];
sx q[3];
rz(-0.20874617) q[3];
sx q[3];
rz(1.1948208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656972) q[0];
sx q[0];
rz(-2.3493451) q[0];
sx q[0];
rz(-2.3727697) q[0];
rz(-2.5104751) q[1];
sx q[1];
rz(-1.9335856) q[1];
sx q[1];
rz(2.9879976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4174059) q[0];
sx q[0];
rz(-1.0225147) q[0];
sx q[0];
rz(2.7566657) q[0];
x q[1];
rz(0.65122143) q[2];
sx q[2];
rz(-1.792728) q[2];
sx q[2];
rz(2.5188308) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92520095) q[1];
sx q[1];
rz(-1.4633302) q[1];
sx q[1];
rz(0.28285427) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0426703) q[3];
sx q[3];
rz(-2.3677845) q[3];
sx q[3];
rz(-2.9681924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.067387335) q[2];
sx q[2];
rz(-1.5978483) q[2];
sx q[2];
rz(0.81573168) q[2];
rz(3.0873114) q[3];
sx q[3];
rz(-1.7547601) q[3];
sx q[3];
rz(1.3666216) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5197007) q[0];
sx q[0];
rz(-2.0760355) q[0];
sx q[0];
rz(-0.52072293) q[0];
rz(-2.0479274) q[1];
sx q[1];
rz(-0.92898527) q[1];
sx q[1];
rz(-1.3826133) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4224515) q[0];
sx q[0];
rz(-2.1364742) q[0];
sx q[0];
rz(-2.9702787) q[0];
rz(-2.6135234) q[2];
sx q[2];
rz(-1.9667528) q[2];
sx q[2];
rz(1.89324) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0212958) q[1];
sx q[1];
rz(-2.3705934) q[1];
sx q[1];
rz(-0.96012562) q[1];
rz(-2.6723271) q[3];
sx q[3];
rz(-1.7466144) q[3];
sx q[3];
rz(-1.8213391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9426721) q[2];
sx q[2];
rz(-2.7896546) q[2];
sx q[2];
rz(1.4658296) q[2];
rz(-0.28907019) q[3];
sx q[3];
rz(-1.7710268) q[3];
sx q[3];
rz(1.2808965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11956231) q[0];
sx q[0];
rz(-0.95198315) q[0];
sx q[0];
rz(-2.6493678) q[0];
rz(0.64552632) q[1];
sx q[1];
rz(-1.5954834) q[1];
sx q[1];
rz(-2.2135977) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4700025) q[0];
sx q[0];
rz(-0.3151463) q[0];
sx q[0];
rz(1.1130087) q[0];
rz(2.7828752) q[2];
sx q[2];
rz(-1.2041098) q[2];
sx q[2];
rz(2.3334954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2454353) q[1];
sx q[1];
rz(-2.3801374) q[1];
sx q[1];
rz(1.1767964) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8230439) q[3];
sx q[3];
rz(-1.2179228) q[3];
sx q[3];
rz(-0.039963756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8202028) q[2];
sx q[2];
rz(-1.6150183) q[2];
sx q[2];
rz(2.3825633) q[2];
rz(-1.2718893) q[3];
sx q[3];
rz(-1.3262409) q[3];
sx q[3];
rz(-0.71182865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3439381) q[0];
sx q[0];
rz(-0.54867083) q[0];
sx q[0];
rz(0.77242533) q[0];
rz(1.7732357) q[1];
sx q[1];
rz(-2.0388347) q[1];
sx q[1];
rz(-0.30430421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876831) q[0];
sx q[0];
rz(-1.4178313) q[0];
sx q[0];
rz(1.3970333) q[0];
rz(-pi) q[1];
rz(-1.7965806) q[2];
sx q[2];
rz(-1.5743557) q[2];
sx q[2];
rz(-0.6706711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58271839) q[1];
sx q[1];
rz(-2.3811023) q[1];
sx q[1];
rz(0.26207102) q[1];
x q[2];
rz(2.7469962) q[3];
sx q[3];
rz(-1.6871042) q[3];
sx q[3];
rz(-1.1726086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.407939) q[2];
sx q[2];
rz(-0.38241461) q[2];
sx q[2];
rz(-1.1747053) q[2];
rz(2.5041653) q[3];
sx q[3];
rz(-1.2624242) q[3];
sx q[3];
rz(-1.9185965) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735073) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(0.77274957) q[0];
rz(-2.6840774) q[1];
sx q[1];
rz(-1.3950149) q[1];
sx q[1];
rz(-1.7083907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.820571) q[0];
sx q[0];
rz(-2.1863156) q[0];
sx q[0];
rz(0.55900283) q[0];
rz(3.0728293) q[2];
sx q[2];
rz(-1.9330852) q[2];
sx q[2];
rz(0.96093528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6683176) q[1];
sx q[1];
rz(-0.5328446) q[1];
sx q[1];
rz(-2.1168296) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97416157) q[3];
sx q[3];
rz(-1.0785558) q[3];
sx q[3];
rz(-0.31271586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92051238) q[2];
sx q[2];
rz(-2.1996193) q[2];
sx q[2];
rz(2.4648049) q[2];
rz(-1.5306728) q[3];
sx q[3];
rz(-2.8407606) q[3];
sx q[3];
rz(1.0845186) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21657011) q[0];
sx q[0];
rz(-1.3022447) q[0];
sx q[0];
rz(-1.362823) q[0];
rz(-2.2577747) q[1];
sx q[1];
rz(-2.4088036) q[1];
sx q[1];
rz(2.1688681) q[1];
rz(2.6893546) q[2];
sx q[2];
rz(-2.0932719) q[2];
sx q[2];
rz(-1.1050638) q[2];
rz(1.3317713) q[3];
sx q[3];
rz(-2.0229133) q[3];
sx q[3];
rz(-2.3537707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
