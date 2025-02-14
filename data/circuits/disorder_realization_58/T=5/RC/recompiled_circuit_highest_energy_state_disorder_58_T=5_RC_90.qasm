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
rz(-2.8623924) q[0];
sx q[0];
rz(-1.8912127) q[0];
sx q[0];
rz(-3.1412636) q[0];
rz(1.1286796) q[1];
sx q[1];
rz(5.1976701) q[1];
sx q[1];
rz(12.357098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34680172) q[0];
sx q[0];
rz(-2.8626094) q[0];
sx q[0];
rz(-2.1643815) q[0];
x q[1];
rz(2.5821788) q[2];
sx q[2];
rz(-1.4296247) q[2];
sx q[2];
rz(0.57746294) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1118374) q[1];
sx q[1];
rz(-1.8833813) q[1];
sx q[1];
rz(-1.8769911) q[1];
x q[2];
rz(-2.9023323) q[3];
sx q[3];
rz(-1.3321028) q[3];
sx q[3];
rz(-2.8857114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2613232) q[2];
sx q[2];
rz(-1.3618733) q[2];
sx q[2];
rz(0.23639354) q[2];
rz(0.38293019) q[3];
sx q[3];
rz(-1.087944) q[3];
sx q[3];
rz(1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3799389) q[0];
sx q[0];
rz(-2.8528657) q[0];
sx q[0];
rz(0.54288236) q[0];
rz(-0.44034964) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(-2.0139258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5701582) q[0];
sx q[0];
rz(-1.2255403) q[0];
sx q[0];
rz(-0.66988129) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2067546) q[2];
sx q[2];
rz(-1.4439772) q[2];
sx q[2];
rz(-0.93709125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.013036916) q[1];
sx q[1];
rz(-0.77289509) q[1];
sx q[1];
rz(1.7835288) q[1];
rz(-pi) q[2];
rz(-1.1912958) q[3];
sx q[3];
rz(-0.89794011) q[3];
sx q[3];
rz(-1.4899474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27309624) q[2];
sx q[2];
rz(-2.345583) q[2];
sx q[2];
rz(1.2924755) q[2];
rz(-2.3097307) q[3];
sx q[3];
rz(-1.4757194) q[3];
sx q[3];
rz(2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6865987) q[0];
sx q[0];
rz(-2.7092777) q[0];
sx q[0];
rz(-1.5469714) q[0];
rz(0.0066241344) q[1];
sx q[1];
rz(-2.6057656) q[1];
sx q[1];
rz(-2.5033902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6168882) q[0];
sx q[0];
rz(-2.0456438) q[0];
sx q[0];
rz(3.1010702) q[0];
x q[1];
rz(1.9784967) q[2];
sx q[2];
rz(-0.8434274) q[2];
sx q[2];
rz(-0.076956017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3824752) q[1];
sx q[1];
rz(-0.51847547) q[1];
sx q[1];
rz(0.7722968) q[1];
rz(-pi) q[2];
rz(2.1177111) q[3];
sx q[3];
rz(-1.7555093) q[3];
sx q[3];
rz(0.46177542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39614761) q[2];
sx q[2];
rz(-3.1148995) q[2];
sx q[2];
rz(-2.2849582) q[2];
rz(-2.4062697) q[3];
sx q[3];
rz(-1.4152799) q[3];
sx q[3];
rz(2.0935238) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51744866) q[0];
sx q[0];
rz(-0.051001661) q[0];
sx q[0];
rz(2.2080102) q[0];
rz(-2.3816662) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(-1.5911721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55465126) q[0];
sx q[0];
rz(-1.4822032) q[0];
sx q[0];
rz(-2.8604126) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43689219) q[2];
sx q[2];
rz(-1.5073188) q[2];
sx q[2];
rz(1.2530099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5129169) q[1];
sx q[1];
rz(-2.4077291) q[1];
sx q[1];
rz(-0.34725125) q[1];
rz(-pi) q[2];
rz(-1.4099259) q[3];
sx q[3];
rz(-1.817559) q[3];
sx q[3];
rz(-2.3575229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4178702) q[2];
sx q[2];
rz(-1.3891862) q[2];
sx q[2];
rz(2.1088481) q[2];
rz(-2.6241153) q[3];
sx q[3];
rz(-1.4314707) q[3];
sx q[3];
rz(-1.8364068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378157) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(-1.284449) q[0];
rz(2.2880554) q[1];
sx q[1];
rz(-1.4505016) q[1];
sx q[1];
rz(0.4761214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.170661) q[0];
sx q[0];
rz(-1.8722152) q[0];
sx q[0];
rz(2.0975251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2849152) q[2];
sx q[2];
rz(-1.21884) q[2];
sx q[2];
rz(-0.98796036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13083982) q[1];
sx q[1];
rz(-0.71341842) q[1];
sx q[1];
rz(2.6793007) q[1];
rz(-pi) q[2];
rz(3.0689803) q[3];
sx q[3];
rz(-2.3092306) q[3];
sx q[3];
rz(1.2439963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7242754) q[2];
sx q[2];
rz(-0.83160916) q[2];
sx q[2];
rz(1.5403436) q[2];
rz(-0.075751461) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(2.4116404) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1002355) q[0];
sx q[0];
rz(-2.5900109) q[0];
sx q[0];
rz(-1.7919737) q[0];
rz(1.8655818) q[1];
sx q[1];
rz(-0.7834692) q[1];
sx q[1];
rz(1.0412201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2305377) q[0];
sx q[0];
rz(-1.7887702) q[0];
sx q[0];
rz(-0.75454069) q[0];
x q[1];
rz(-3.1026672) q[2];
sx q[2];
rz(-2.3934529) q[2];
sx q[2];
rz(0.85252658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0184082) q[1];
sx q[1];
rz(-0.45547541) q[1];
sx q[1];
rz(0.60103086) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3833748) q[3];
sx q[3];
rz(-1.3617523) q[3];
sx q[3];
rz(-0.93702873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87018806) q[2];
sx q[2];
rz(-1.9726334) q[2];
sx q[2];
rz(0.99883396) q[2];
rz(0.53423229) q[3];
sx q[3];
rz(-1.9637008) q[3];
sx q[3];
rz(-1.4786725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6967195) q[0];
sx q[0];
rz(-0.98170009) q[0];
sx q[0];
rz(1.556742) q[0];
rz(-2.059767) q[1];
sx q[1];
rz(-1.6422681) q[1];
sx q[1];
rz(-0.36344847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7963821) q[0];
sx q[0];
rz(-1.78336) q[0];
sx q[0];
rz(1.5770771) q[0];
x q[1];
rz(0.3360904) q[2];
sx q[2];
rz(-2.828183) q[2];
sx q[2];
rz(-1.4614507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9103437) q[1];
sx q[1];
rz(-1.1414613) q[1];
sx q[1];
rz(2.3554159) q[1];
x q[2];
rz(-2.0601252) q[3];
sx q[3];
rz(-0.31144588) q[3];
sx q[3];
rz(1.2490727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5536559) q[2];
sx q[2];
rz(-2.9755972) q[2];
sx q[2];
rz(-2.0641649) q[2];
rz(1.59498) q[3];
sx q[3];
rz(-1.8190106) q[3];
sx q[3];
rz(0.67302978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308291) q[0];
sx q[0];
rz(-1.4973233) q[0];
sx q[0];
rz(2.8610435) q[0];
rz(-2.5827017) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(1.9299318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9292363) q[0];
sx q[0];
rz(-0.87221891) q[0];
sx q[0];
rz(0.86232604) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40139426) q[2];
sx q[2];
rz(-0.52782431) q[2];
sx q[2];
rz(-2.2032705) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9338464) q[1];
sx q[1];
rz(-0.7911111) q[1];
sx q[1];
rz(-1.6027811) q[1];
x q[2];
rz(-1.5455963) q[3];
sx q[3];
rz(-0.19807252) q[3];
sx q[3];
rz(0.36626178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1739379) q[2];
sx q[2];
rz(-1.329198) q[2];
sx q[2];
rz(-2.3363414) q[2];
rz(1.1236745) q[3];
sx q[3];
rz(-2.0011963) q[3];
sx q[3];
rz(-2.8203188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6174018) q[0];
sx q[0];
rz(-0.62365714) q[0];
sx q[0];
rz(-0.50530857) q[0];
rz(2.4303923) q[1];
sx q[1];
rz(-0.86831793) q[1];
sx q[1];
rz(-2.0288846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8998535) q[0];
sx q[0];
rz(-2.02075) q[0];
sx q[0];
rz(-1.8983272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10782055) q[2];
sx q[2];
rz(-1.6514213) q[2];
sx q[2];
rz(0.18641678) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1755274) q[1];
sx q[1];
rz(-1.6972145) q[1];
sx q[1];
rz(-1.1299707) q[1];
rz(1.6203513) q[3];
sx q[3];
rz(-0.55618868) q[3];
sx q[3];
rz(-0.24133989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6342371) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(-0.235454) q[2];
rz(-0.37774751) q[3];
sx q[3];
rz(-2.7538959) q[3];
sx q[3];
rz(-0.050067576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560051) q[0];
sx q[0];
rz(-1.2879141) q[0];
sx q[0];
rz(-0.054314286) q[0];
rz(2.3771225) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(2.4321709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55290907) q[0];
sx q[0];
rz(-1.3257265) q[0];
sx q[0];
rz(-3.0348196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26086678) q[2];
sx q[2];
rz(-1.7013019) q[2];
sx q[2];
rz(2.5640287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1349984) q[1];
sx q[1];
rz(-1.6853764) q[1];
sx q[1];
rz(2.5496818) q[1];
rz(-1.4414431) q[3];
sx q[3];
rz(-1.23063) q[3];
sx q[3];
rz(0.29290712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28598049) q[2];
sx q[2];
rz(-1.4748272) q[2];
sx q[2];
rz(-3.0484071) q[2];
rz(-0.5992254) q[3];
sx q[3];
rz(-2.5803876) q[3];
sx q[3];
rz(2.3563201) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7080606) q[0];
sx q[0];
rz(-2.1716433) q[0];
sx q[0];
rz(-0.94950983) q[0];
rz(0.064619725) q[1];
sx q[1];
rz(-0.043793543) q[1];
sx q[1];
rz(0.66088062) q[1];
rz(1.0192766) q[2];
sx q[2];
rz(-1.6342155) q[2];
sx q[2];
rz(-1.8104959) q[2];
rz(-1.384769) q[3];
sx q[3];
rz(-0.91384619) q[3];
sx q[3];
rz(-2.1887907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
