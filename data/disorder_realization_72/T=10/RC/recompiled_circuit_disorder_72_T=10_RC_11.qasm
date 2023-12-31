OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(-1.6367903) q[0];
rz(-pi) q[1];
rz(0.84988611) q[2];
sx q[2];
rz(-1.4404785) q[2];
sx q[2];
rz(2.7352114) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5159113) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(0.73966571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-0.10676521) q[3];
sx q[3];
rz(-1.8330542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(-2.1271558) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-1.0043253) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(2.1388334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28264499) q[0];
sx q[0];
rz(-0.98416057) q[0];
sx q[0];
rz(-1.0413917) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3833212) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(-1.0152917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.037307449) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(2.9557455) q[1];
rz(-pi) q[2];
rz(-0.49591222) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944091) q[0];
sx q[0];
rz(-0.43617019) q[0];
sx q[0];
rz(2.7437074) q[0];
rz(-pi) q[1];
rz(-1.1728889) q[2];
sx q[2];
rz(-1.8430317) q[2];
sx q[2];
rz(0.18323252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5933697) q[1];
sx q[1];
rz(1.2858461) q[1];
x q[2];
rz(-1.4530573) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(-1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94565369) q[0];
sx q[0];
rz(-1.6473856) q[0];
sx q[0];
rz(-0.73777317) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29454622) q[2];
sx q[2];
rz(-1.4809594) q[2];
sx q[2];
rz(-1.3603269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6039227) q[1];
sx q[1];
rz(-2.1426755) q[1];
sx q[1];
rz(-0.13725431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6226193) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(0.94829544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0147858) q[0];
sx q[0];
rz(-0.94313398) q[0];
sx q[0];
rz(0.058237596) q[0];
x q[1];
rz(1.615633) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(-2.5068138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.535714) q[1];
sx q[1];
rz(-1.6026346) q[1];
sx q[1];
rz(1.72623) q[1];
rz(-pi) q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-0.92901232) q[3];
sx q[3];
rz(-1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(0.92818964) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2561589) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(-2.8594349) q[0];
x q[1];
rz(1.2601389) q[2];
sx q[2];
rz(-2.0271745) q[2];
sx q[2];
rz(-3.0820456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0551758) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(3.0228826) q[1];
rz(-pi) q[2];
rz(-2.5572705) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(0.82715183) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-1.126359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(0.29445946) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14234219) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(3.1359283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4451051) q[1];
sx q[1];
rz(-1.3828053) q[1];
sx q[1];
rz(0.45447116) q[1];
rz(-pi) q[2];
rz(1.2231636) q[3];
sx q[3];
rz(-1.6602483) q[3];
sx q[3];
rz(1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(-2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643124) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(2.2420922) q[0];
x q[1];
rz(-1.007349) q[2];
sx q[2];
rz(-1.9512366) q[2];
sx q[2];
rz(0.13605875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9084839) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(-1.7912346) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(2.5194871) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(0.73293066) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(2.656235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982614) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(-1.5166548) q[0];
x q[1];
rz(1.4114686) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(-1.9221905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-0.81307462) q[1];
x q[2];
rz(2.5599307) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(0.5982582) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(-1.422545) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(-1.460618) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.9365063) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7353021) q[3];
sx q[3];
rz(-1.7269772) q[3];
sx q[3];
rz(-2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-2.5219953) q[2];
sx q[2];
rz(-1.6710812) q[2];
sx q[2];
rz(-2.8166213) q[2];
rz(2.2451154) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
