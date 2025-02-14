OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(1.6311128) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(4.4463867) q[1];
sx q[1];
rz(10.211853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2505456) q[0];
sx q[0];
rz(-0.68792397) q[0];
sx q[0];
rz(2.964818) q[0];
rz(0.39648367) q[2];
sx q[2];
rz(-1.6210437) q[2];
sx q[2];
rz(-0.32279245) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30336994) q[1];
sx q[1];
rz(-1.359424) q[1];
sx q[1];
rz(-1.6081393) q[1];
rz(-pi) q[2];
rz(2.1605499) q[3];
sx q[3];
rz(-2.1043805) q[3];
sx q[3];
rz(-0.55280815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2934908) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(-0.21609406) q[2];
rz(-2.6381093) q[3];
sx q[3];
rz(-1.8899906) q[3];
sx q[3];
rz(3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(2.766372) q[0];
rz(-0.61620617) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(0.18888758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096029671) q[0];
sx q[0];
rz(-2.3052373) q[0];
sx q[0];
rz(2.1406853) q[0];
rz(-pi) q[1];
rz(2.1470931) q[2];
sx q[2];
rz(-2.3220571) q[2];
sx q[2];
rz(-0.34162765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0068854) q[1];
sx q[1];
rz(-1.3682248) q[1];
sx q[1];
rz(2.2565875) q[1];
x q[2];
rz(2.5152605) q[3];
sx q[3];
rz(-0.62343684) q[3];
sx q[3];
rz(-1.745519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6836493) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(-1.1434435) q[2];
rz(-2.5055366) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77883333) q[0];
sx q[0];
rz(-0.23985671) q[0];
sx q[0];
rz(1.1676316) q[0];
rz(3.0150343) q[1];
sx q[1];
rz(-1.626222) q[1];
sx q[1];
rz(-2.5680241) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48808778) q[0];
sx q[0];
rz(-1.4801985) q[0];
sx q[0];
rz(0.82838627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5218256) q[2];
sx q[2];
rz(-2.7006442) q[2];
sx q[2];
rz(2.3201705) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14309569) q[1];
sx q[1];
rz(-1.6187889) q[1];
sx q[1];
rz(2.001862) q[1];
x q[2];
rz(-2.6010547) q[3];
sx q[3];
rz(-0.73523587) q[3];
sx q[3];
rz(-1.0395887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6557189) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(1.5029079) q[2];
rz(1.2949299) q[3];
sx q[3];
rz(-2.86627) q[3];
sx q[3];
rz(2.3143229) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5591705) q[0];
sx q[0];
rz(-2.9952413) q[0];
sx q[0];
rz(2.225112) q[0];
rz(0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(-0.697335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50619102) q[0];
sx q[0];
rz(-0.40987337) q[0];
sx q[0];
rz(-0.91110595) q[0];
rz(-pi) q[1];
rz(-0.38390145) q[2];
sx q[2];
rz(-2.4426201) q[2];
sx q[2];
rz(-1.1363824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26665822) q[1];
sx q[1];
rz(-1.8968079) q[1];
sx q[1];
rz(-0.74649109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.55682) q[3];
sx q[3];
rz(-1.7009354) q[3];
sx q[3];
rz(1.1443477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9852898) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(2.018833) q[2];
rz(-1.5924234) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(2.7482225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(1.1291809) q[0];
rz(-0.87942266) q[1];
sx q[1];
rz(-1.3112473) q[1];
sx q[1];
rz(0.62854356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56282069) q[0];
sx q[0];
rz(-1.6878753) q[0];
sx q[0];
rz(-1.4463918) q[0];
rz(-pi) q[1];
rz(-2.4714575) q[2];
sx q[2];
rz(-2.3628855) q[2];
sx q[2];
rz(0.35691945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6595093) q[1];
sx q[1];
rz(-2.1225559) q[1];
sx q[1];
rz(0.047288069) q[1];
rz(-pi) q[2];
rz(0.83129779) q[3];
sx q[3];
rz(-1.4834838) q[3];
sx q[3];
rz(-2.1082474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(1.3735695) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-1.0765358) q[3];
sx q[3];
rz(-0.086418644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81109989) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(-2.9789441) q[0];
rz(1.9074408) q[1];
sx q[1];
rz(-1.994543) q[1];
sx q[1];
rz(-0.92076921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5266357) q[0];
sx q[0];
rz(-1.6713854) q[0];
sx q[0];
rz(1.3822162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55947907) q[2];
sx q[2];
rz(-2.1396416) q[2];
sx q[2];
rz(0.43348962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28345767) q[1];
sx q[1];
rz(-0.75443903) q[1];
sx q[1];
rz(0.094875022) q[1];
rz(-pi) q[2];
rz(0.1699888) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(1.9758177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84772253) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(-3.0976307) q[2];
rz(-1.4219159) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(-0.034984263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1304355) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-3.0479808) q[0];
rz(-1.0253819) q[1];
sx q[1];
rz(-1.5954115) q[1];
sx q[1];
rz(-2.3536918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9196531) q[0];
sx q[0];
rz(-1.0353695) q[0];
sx q[0];
rz(0.8485283) q[0];
rz(-pi) q[1];
rz(-1.4008994) q[2];
sx q[2];
rz(-2.133103) q[2];
sx q[2];
rz(0.0058836939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8389143) q[1];
sx q[1];
rz(-0.60015772) q[1];
sx q[1];
rz(-3.0693574) q[1];
rz(-pi) q[2];
rz(2.2719194) q[3];
sx q[3];
rz(-0.18571412) q[3];
sx q[3];
rz(2.7839422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0755783) q[2];
sx q[2];
rz(-2.5983073) q[2];
sx q[2];
rz(2.1755966) q[2];
rz(0.14723369) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(-2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(1.8889486) q[0];
rz(-0.057627536) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(2.7431814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1712043) q[0];
sx q[0];
rz(-0.96231595) q[0];
sx q[0];
rz(1.7826016) q[0];
rz(-pi) q[1];
rz(1.4008825) q[2];
sx q[2];
rz(-2.0644098) q[2];
sx q[2];
rz(2.5371671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9908473) q[1];
sx q[1];
rz(-1.3708893) q[1];
sx q[1];
rz(-3.1374817) q[1];
rz(-pi) q[2];
rz(-0.54315217) q[3];
sx q[3];
rz(-1.6622322) q[3];
sx q[3];
rz(0.1426129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5742699) q[2];
sx q[2];
rz(-1.2719354) q[2];
sx q[2];
rz(2.2197913) q[2];
rz(-2.4343893) q[3];
sx q[3];
rz(-1.0450398) q[3];
sx q[3];
rz(-0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(0.90712547) q[0];
rz(2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(-0.60636806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7146695) q[0];
sx q[0];
rz(-1.4703317) q[0];
sx q[0];
rz(0.31025525) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81524088) q[2];
sx q[2];
rz(-0.70910197) q[2];
sx q[2];
rz(0.39295024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38708147) q[1];
sx q[1];
rz(-2.3007327) q[1];
sx q[1];
rz(0.39772268) q[1];
rz(-pi) q[2];
rz(1.4616897) q[3];
sx q[3];
rz(-2.2266995) q[3];
sx q[3];
rz(0.65306585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85764641) q[2];
sx q[2];
rz(-1.593677) q[2];
sx q[2];
rz(0.66961163) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-2.0958869) q[3];
sx q[3];
rz(1.7414352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46212101) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(2.2456428) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(-1.5641854) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58916726) q[0];
sx q[0];
rz(-1.5755113) q[0];
sx q[0];
rz(1.5472359) q[0];
rz(0.51367398) q[2];
sx q[2];
rz(-2.3875055) q[2];
sx q[2];
rz(-2.9083792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60301723) q[1];
sx q[1];
rz(-1.9780206) q[1];
sx q[1];
rz(-1.0124221) q[1];
rz(-1.0684929) q[3];
sx q[3];
rz(-2.4728259) q[3];
sx q[3];
rz(-2.7947102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3080052) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(1.8316815) q[2];
rz(-1.1854019) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014521) q[0];
sx q[0];
rz(-1.3012713) q[0];
sx q[0];
rz(-0.95300994) q[0];
rz(-0.078527191) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-2.2453515) q[2];
sx q[2];
rz(-1.7755388) q[2];
sx q[2];
rz(-0.15646738) q[2];
rz(2.9361712) q[3];
sx q[3];
rz(-1.4170437) q[3];
sx q[3];
rz(0.99544256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
