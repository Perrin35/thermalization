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
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5242033) q[0];
sx q[0];
rz(-0.19771068) q[0];
sx q[0];
rz(0.80394216) q[0];
x q[1];
rz(0.013597537) q[2];
sx q[2];
rz(-2.6981079) q[2];
sx q[2];
rz(0.19946239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9182723) q[1];
sx q[1];
rz(-1.2967355) q[1];
sx q[1];
rz(2.5943701) q[1];
rz(-pi) q[2];
rz(-0.080022607) q[3];
sx q[3];
rz(-1.140174) q[3];
sx q[3];
rz(2.9324556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67285377) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(-1.2057745) q[2];
rz(-2.2852211) q[3];
sx q[3];
rz(-0.93256336) q[3];
sx q[3];
rz(-0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641651) q[0];
sx q[0];
rz(-0.39746118) q[0];
sx q[0];
rz(-0.20832668) q[0];
rz(-1.2267998) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(2.5770381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3608332) q[0];
sx q[0];
rz(-1.5281023) q[0];
sx q[0];
rz(1.4111931) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54148285) q[2];
sx q[2];
rz(-2.3332806) q[2];
sx q[2];
rz(1.9549416) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7617088) q[1];
sx q[1];
rz(-1.392388) q[1];
sx q[1];
rz(0.70990035) q[1];
x q[2];
rz(1.0309593) q[3];
sx q[3];
rz(-1.5542277) q[3];
sx q[3];
rz(1.4318717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.2496639) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(-2.6460323) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(-1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946852) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(1.2184628) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(2.7168435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0828511) q[0];
sx q[0];
rz(-0.070151873) q[0];
sx q[0];
rz(-2.6741316) q[0];
rz(-pi) q[1];
rz(0.23992691) q[2];
sx q[2];
rz(-1.3245965) q[2];
sx q[2];
rz(2.4407516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9806165) q[1];
sx q[1];
rz(-0.38432877) q[1];
sx q[1];
rz(2.9381782) q[1];
rz(-pi) q[2];
rz(3.0851641) q[3];
sx q[3];
rz(-0.61739576) q[3];
sx q[3];
rz(1.4790725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4212627) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(3.1171411) q[2];
rz(1.9931741) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-0.20632437) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(0.21388737) q[0];
rz(-2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-2.4549386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7711597) q[0];
sx q[0];
rz(-0.70777445) q[0];
sx q[0];
rz(0.64162713) q[0];
x q[1];
rz(-1.114421) q[2];
sx q[2];
rz(-0.59774146) q[2];
sx q[2];
rz(-0.98243143) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(-0.13607009) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88249256) q[3];
sx q[3];
rz(-0.63418856) q[3];
sx q[3];
rz(-0.89804441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.116918) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(2.4521258) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(-0.57233468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6784742) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(-1.0372739) q[0];
rz(0.67266881) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(-0.79576463) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507495) q[0];
sx q[0];
rz(-1.8863057) q[0];
sx q[0];
rz(-2.1054224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8383083) q[2];
sx q[2];
rz(-2.1223084) q[2];
sx q[2];
rz(1.9818652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63043609) q[1];
sx q[1];
rz(-1.9025584) q[1];
sx q[1];
rz(-2.3219882) q[1];
rz(-0.66821258) q[3];
sx q[3];
rz(-1.8132121) q[3];
sx q[3];
rz(-0.53430623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(-2.3058092) q[2];
rz(-0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68818727) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(-1.2123464) q[0];
rz(1.0351828) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(-1.0119247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31549227) q[0];
sx q[0];
rz(-1.4619215) q[0];
sx q[0];
rz(1.3988162) q[0];
x q[1];
rz(3.1306562) q[2];
sx q[2];
rz(-1.2669529) q[2];
sx q[2];
rz(2.2774709) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3287836) q[1];
sx q[1];
rz(-1.145383) q[1];
sx q[1];
rz(1.119647) q[1];
rz(-pi) q[2];
rz(2.7583241) q[3];
sx q[3];
rz(-1.82194) q[3];
sx q[3];
rz(-2.0519837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(2.193023) q[2];
rz(-2.6712724) q[3];
sx q[3];
rz(-1.0985273) q[3];
sx q[3];
rz(-1.8868014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5941641) q[0];
sx q[0];
rz(-1.0110039) q[0];
sx q[0];
rz(-1.2642566) q[0];
rz(-2.3612379) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(2.9497214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.980968) q[0];
sx q[0];
rz(-1.1066529) q[0];
sx q[0];
rz(-1.1726517) q[0];
rz(-pi) q[1];
rz(2.3290995) q[2];
sx q[2];
rz(-2.6408962) q[2];
sx q[2];
rz(2.814584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7765155) q[1];
sx q[1];
rz(-2.0151225) q[1];
sx q[1];
rz(-1.278484) q[1];
rz(-pi) q[2];
x q[2];
rz(2.170788) q[3];
sx q[3];
rz(-1.8626584) q[3];
sx q[3];
rz(1.4960904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(-0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845487) q[0];
sx q[0];
rz(-0.024024155) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(1.7441033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5834806) q[0];
sx q[0];
rz(-1.8761823) q[0];
sx q[0];
rz(-2.7598513) q[0];
rz(-pi) q[1];
rz(-1.0447803) q[2];
sx q[2];
rz(-0.57990852) q[2];
sx q[2];
rz(-1.2666262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8694526) q[1];
sx q[1];
rz(-1.6873797) q[1];
sx q[1];
rz(-3.0805853) q[1];
rz(-2.0040183) q[3];
sx q[3];
rz(-0.94846359) q[3];
sx q[3];
rz(-2.4985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(-1.5020471) q[2];
rz(-1.8223193) q[3];
sx q[3];
rz(-0.84984142) q[3];
sx q[3];
rz(1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424292) q[0];
sx q[0];
rz(-2.2971239) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(-0.49081048) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(-1.6455654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67453066) q[0];
sx q[0];
rz(-0.26755324) q[0];
sx q[0];
rz(-0.48488276) q[0];
rz(1.7424351) q[2];
sx q[2];
rz(-1.9783859) q[2];
sx q[2];
rz(1.7625347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89801187) q[1];
sx q[1];
rz(-1.3286546) q[1];
sx q[1];
rz(0.33278521) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4946901) q[3];
sx q[3];
rz(-2.1015014) q[3];
sx q[3];
rz(-0.87065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8533123) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(1.3221928) q[2];
rz(2.7754122) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.2243097) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(-3.068058) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(1.2988466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.038485) q[0];
sx q[0];
rz(-1.847205) q[0];
sx q[0];
rz(1.7544708) q[0];
rz(-pi) q[1];
rz(-1.0825363) q[2];
sx q[2];
rz(-1.6886782) q[2];
sx q[2];
rz(-2.0669017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9292752) q[1];
sx q[1];
rz(-1.8061418) q[1];
sx q[1];
rz(1.9566243) q[1];
rz(-0.11839827) q[3];
sx q[3];
rz(-2.3984635) q[3];
sx q[3];
rz(1.2606674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1765882) q[2];
sx q[2];
rz(-0.95866385) q[2];
sx q[2];
rz(-2.5625572) q[2];
rz(1.6252919) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(-1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2753006) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(0.91118377) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(-1.5742792) q[2];
sx q[2];
rz(-1.9949253) q[2];
sx q[2];
rz(-0.16535769) q[2];
rz(-0.95355036) q[3];
sx q[3];
rz(-1.8051157) q[3];
sx q[3];
rz(2.3122862) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
