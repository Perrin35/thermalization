OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53579522) q[0];
sx q[0];
rz(-1.3942766) q[0];
sx q[0];
rz(-1.4730886) q[0];
rz(-pi) q[1];
rz(-1.2486357) q[2];
sx q[2];
rz(-0.82167168) q[2];
sx q[2];
rz(1.2920213) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40542422) q[1];
sx q[1];
rz(-2.5129457) q[1];
sx q[1];
rz(2.5581762) q[1];
rz(1.3084859) q[3];
sx q[3];
rz(-2.9648818) q[3];
sx q[3];
rz(1.7929329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9464232) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(0.75631022) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-0.38683495) q[0];
rz(2.6392) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5997255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8097336) q[0];
sx q[0];
rz(-0.74685687) q[0];
sx q[0];
rz(-2.5144308) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1142271) q[2];
sx q[2];
rz(-1.7184988) q[2];
sx q[2];
rz(2.4545124) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7911885) q[1];
sx q[1];
rz(-1.5608619) q[1];
sx q[1];
rz(-3.105858) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76962556) q[3];
sx q[3];
rz(-2.5105021) q[3];
sx q[3];
rz(-2.1345994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(-1.4146457) q[2];
rz(2.291262) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.458228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787815) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(-0.61022726) q[0];
rz(-1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(2.1496225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45527601) q[0];
sx q[0];
rz(-1.0193829) q[0];
sx q[0];
rz(-0.21689143) q[0];
x q[1];
rz(1.3685554) q[2];
sx q[2];
rz(-2.2291406) q[2];
sx q[2];
rz(-2.1683247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75533463) q[1];
sx q[1];
rz(-1.4812246) q[1];
sx q[1];
rz(-2.7002525) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3258341) q[3];
sx q[3];
rz(-1.1253329) q[3];
sx q[3];
rz(2.9374292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38301864) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(-2.7475083) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(2.7365141) q[0];
rz(-0.69008094) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(-2.4437723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9245802) q[0];
sx q[0];
rz(-3.0487061) q[0];
sx q[0];
rz(-1.8266982) q[0];
rz(0.54450808) q[2];
sx q[2];
rz(-1.240057) q[2];
sx q[2];
rz(-0.43560189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3623558) q[1];
sx q[1];
rz(-2.6168565) q[1];
sx q[1];
rz(-1.1974105) q[1];
x q[2];
rz(1.9911733) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(0.099893173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.7439338) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(-2.737282) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25800911) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(-2.1160545) q[0];
rz(2.569596) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(-0.62932032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176661) q[0];
sx q[0];
rz(-1.1610982) q[0];
sx q[0];
rz(1.2439338) q[0];
rz(-pi) q[1];
rz(0.25482486) q[2];
sx q[2];
rz(-0.91081589) q[2];
sx q[2];
rz(0.13435907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8992936) q[1];
sx q[1];
rz(-2.0819903) q[1];
sx q[1];
rz(-1.7023229) q[1];
x q[2];
rz(0.96962813) q[3];
sx q[3];
rz(-1.3519577) q[3];
sx q[3];
rz(-0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59262529) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(-1.3458378) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(-0.19601823) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(-2.956399) q[0];
rz(1.7350896) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.3669744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0104116) q[0];
sx q[0];
rz(-0.62958065) q[0];
sx q[0];
rz(-2.7601943) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84080266) q[2];
sx q[2];
rz(-1.3422988) q[2];
sx q[2];
rz(2.218354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1320912) q[1];
sx q[1];
rz(-2.9638634) q[1];
sx q[1];
rz(1.0264261) q[1];
rz(-0.56914764) q[3];
sx q[3];
rz(-1.9372809) q[3];
sx q[3];
rz(1.6397427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(-0.883376) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(-0.81378716) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25154034) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(1.863377) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.7657123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517075) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(2.7566107) q[0];
rz(-pi) q[1];
rz(-1.4543578) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(1.4637228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5232877) q[1];
sx q[1];
rz(-1.4555281) q[1];
sx q[1];
rz(1.4646039) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0828956) q[3];
sx q[3];
rz(-1.241063) q[3];
sx q[3];
rz(-1.2122455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4574796) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(0.11080065) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(-0.69537648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(2.5336174) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(2.9072993) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5406571) q[0];
sx q[0];
rz(-1.0957076) q[0];
sx q[0];
rz(-0.29151543) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4085326) q[2];
sx q[2];
rz(-0.86793938) q[2];
sx q[2];
rz(1.5066063) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5768847) q[1];
sx q[1];
rz(-1.7600093) q[1];
sx q[1];
rz(2.3343759) q[1];
x q[2];
rz(2.0571363) q[3];
sx q[3];
rz(-2.497695) q[3];
sx q[3];
rz(2.7542496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12604788) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(1.4253433) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.4956168) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(-2.1967922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1856857) q[0];
sx q[0];
rz(-2.3392896) q[0];
sx q[0];
rz(0.54922744) q[0];
rz(1.1759042) q[2];
sx q[2];
rz(-1.8443622) q[2];
sx q[2];
rz(-3.0964031) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8560801) q[1];
sx q[1];
rz(-1.3998919) q[1];
sx q[1];
rz(-2.8588365) q[1];
x q[2];
rz(-1.9267843) q[3];
sx q[3];
rz(-1.755852) q[3];
sx q[3];
rz(1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21645674) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(-2.2286041) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(2.250681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83715314) q[0];
sx q[0];
rz(-1.3225978) q[0];
sx q[0];
rz(-2.7960143) q[0];
x q[1];
rz(-2.183379) q[2];
sx q[2];
rz(-1.805086) q[2];
sx q[2];
rz(0.9226375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18894698) q[1];
sx q[1];
rz(-1.9083605) q[1];
sx q[1];
rz(-1.360421) q[1];
rz(-pi) q[2];
rz(1.4952412) q[3];
sx q[3];
rz(-1.75845) q[3];
sx q[3];
rz(-2.5726266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(-1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(0.44617173) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5065153) q[0];
sx q[0];
rz(-1.7000533) q[0];
sx q[0];
rz(0.62361367) q[0];
rz(-2.0093909) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(1.7239465) q[2];
sx q[2];
rz(-1.3071123) q[2];
sx q[2];
rz(1.2109962) q[2];
rz(1.6518456) q[3];
sx q[3];
rz(-2.5935018) q[3];
sx q[3];
rz(0.58542585) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
