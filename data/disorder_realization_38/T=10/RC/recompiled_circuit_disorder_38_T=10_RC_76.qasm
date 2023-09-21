OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2599517) q[0];
sx q[0];
rz(-1.5978442) q[0];
sx q[0];
rz(1.7904439) q[0];
rz(-pi) q[1];
rz(-0.89994853) q[2];
sx q[2];
rz(-2.4953825) q[2];
sx q[2];
rz(-0.83713573) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7632335) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(0.50904973) q[1];
x q[2];
rz(-0.083505587) q[3];
sx q[3];
rz(-2.180035) q[3];
sx q[3];
rz(0.0041675605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(-2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(-2.8711328) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(1.6289904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0892031) q[0];
sx q[0];
rz(-0.69141483) q[0];
sx q[0];
rz(1.0538488) q[0];
x q[1];
rz(-1.5674044) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(1.1018745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43803793) q[1];
sx q[1];
rz(-2.2433271) q[1];
sx q[1];
rz(0.1708252) q[1];
rz(-pi) q[2];
rz(2.2643331) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(-3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495372) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(-1.3489086) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4065811) q[2];
sx q[2];
rz(-1.0501554) q[2];
sx q[2];
rz(-0.28902136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44887603) q[1];
sx q[1];
rz(-1.2926896) q[1];
sx q[1];
rz(2.6443308) q[1];
rz(-pi) q[2];
rz(-0.45332076) q[3];
sx q[3];
rz(-1.3461777) q[3];
sx q[3];
rz(0.9325222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(-2.8682958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647588) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.2980952) q[0];
rz(-0.69584537) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30444333) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(1.3346346) q[1];
rz(-pi) q[2];
rz(0.64951879) q[3];
sx q[3];
rz(-1.4276541) q[3];
sx q[3];
rz(0.79917819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(3.058847) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(-0.98714978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(2.2228918) q[0];
rz(-pi) q[1];
rz(-3.1290595) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(0.24521337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3293216) q[1];
sx q[1];
rz(-1.8835888) q[1];
sx q[1];
rz(2.4815464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9636642) q[3];
sx q[3];
rz(-1.53189) q[3];
sx q[3];
rz(1.9143357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.2794367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50514454) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.5178773) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-2.1589303) q[2];
sx q[2];
rz(2.9786125) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93950677) q[1];
sx q[1];
rz(-2.0923951) q[1];
sx q[1];
rz(3.107198) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1412973) q[3];
sx q[3];
rz(-1.391489) q[3];
sx q[3];
rz(0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-0.77077579) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188748) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-3.0859257) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-3.1380222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(-2.4128782) q[0];
x q[1];
rz(2.3891719) q[2];
sx q[2];
rz(-1.2781029) q[2];
sx q[2];
rz(2.9257286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6047302) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(2.0525949) q[1];
x q[2];
rz(-1.9332063) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(-0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-2.2214831) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(1.2922657) q[0];
rz(0.46531123) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(-1.0516143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.855367) q[1];
sx q[1];
rz(-1.507842) q[1];
sx q[1];
rz(1.0463868) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.802556) q[3];
sx q[3];
rz(-0.44796523) q[3];
sx q[3];
rz(-0.26089222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(2.881799) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(-2.7816539) q[0];
x q[1];
rz(-2.8723268) q[2];
sx q[2];
rz(-0.93062799) q[2];
sx q[2];
rz(-0.38538853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.867013) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(0.32151476) q[1];
rz(-pi) q[2];
rz(2.6960877) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(-1.477369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-0.3607761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4496778) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(1.7382394) q[0];
x q[1];
rz(3.0214494) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-0.10211589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68041486) q[1];
sx q[1];
rz(-0.83592452) q[1];
sx q[1];
rz(1.3780891) q[1];
rz(-pi) q[2];
rz(-2.9721857) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(1.3225079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(2.0104682) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(-0.75136649) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
