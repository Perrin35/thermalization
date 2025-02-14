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
rz(0.063989446) q[0];
sx q[0];
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(1.9995243) q[1];
sx q[1];
rz(-2.6231782) q[1];
sx q[1];
rz(1.4919182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5357852) q[0];
sx q[0];
rz(-1.5174179) q[0];
sx q[0];
rz(-0.018529539) q[0];
rz(-pi) q[1];
rz(2.1989735) q[2];
sx q[2];
rz(-2.16428) q[2];
sx q[2];
rz(-1.3243937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0483889) q[1];
sx q[1];
rz(-1.2902765) q[1];
sx q[1];
rz(-2.2040221) q[1];
x q[2];
rz(0.033871058) q[3];
sx q[3];
rz(-1.3283786) q[3];
sx q[3];
rz(-0.1823547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3492744) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(0.054923687) q[2];
rz(-1.4548291) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(-1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(-0.24800214) q[0];
rz(1.8964881) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.9128333) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354483) q[0];
sx q[0];
rz(-1.7559253) q[0];
sx q[0];
rz(1.6044751) q[0];
rz(-pi) q[1];
rz(2.6284559) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(-3.0233011) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1540888) q[1];
sx q[1];
rz(-2.262907) q[1];
sx q[1];
rz(-1.8471884) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13350962) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(0.51160073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(-0.42932388) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(-2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(-1.6744457) q[0];
rz(-1.4779429) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(2.0416226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29102688) q[0];
sx q[0];
rz(-1.675559) q[0];
sx q[0];
rz(2.9423703) q[0];
rz(2.8520865) q[2];
sx q[2];
rz(-1.2443674) q[2];
sx q[2];
rz(-3.141186) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63640672) q[1];
sx q[1];
rz(-0.51915324) q[1];
sx q[1];
rz(0.7208419) q[1];
rz(0.93020029) q[3];
sx q[3];
rz(-2.4630694) q[3];
sx q[3];
rz(0.92032209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(-1.5477017) q[2];
rz(1.8111546) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(2.0681341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610157) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-2.9837578) q[0];
rz(0.24179587) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(-2.6604624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8932228) q[0];
sx q[0];
rz(-0.71740323) q[0];
sx q[0];
rz(0.1347085) q[0];
rz(-1.3161737) q[2];
sx q[2];
rz(-1.919165) q[2];
sx q[2];
rz(0.35149945) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1930052) q[1];
sx q[1];
rz(-0.95110106) q[1];
sx q[1];
rz(0.42264688) q[1];
x q[2];
rz(-1.8577544) q[3];
sx q[3];
rz(-1.3697512) q[3];
sx q[3];
rz(1.6359117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2699997) q[2];
sx q[2];
rz(-0.81120482) q[2];
sx q[2];
rz(0.2743741) q[2];
rz(-1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6896553) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(-1.0614606) q[0];
rz(2.6929216) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(-1.4400858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3685303) q[0];
sx q[0];
rz(-2.321918) q[0];
sx q[0];
rz(-1.7858265) q[0];
x q[1];
rz(-3.0147047) q[2];
sx q[2];
rz(-1.972162) q[2];
sx q[2];
rz(-0.021406476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60665874) q[1];
sx q[1];
rz(-0.25940653) q[1];
sx q[1];
rz(1.4662798) q[1];
x q[2];
rz(1.7251882) q[3];
sx q[3];
rz(-2.5254211) q[3];
sx q[3];
rz(2.7026724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3638641) q[2];
sx q[2];
rz(-1.8455467) q[2];
sx q[2];
rz(1.9484005) q[2];
rz(2.7423972) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(-0.68049085) q[0];
rz(-2.3091799) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(0.15636538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2966317) q[0];
sx q[0];
rz(-2.8609242) q[0];
sx q[0];
rz(-0.190221) q[0];
x q[1];
rz(-0.6758718) q[2];
sx q[2];
rz(-2.1891433) q[2];
sx q[2];
rz(-2.4175274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.865327) q[1];
sx q[1];
rz(-0.4368096) q[1];
sx q[1];
rz(-2.275009) q[1];
rz(-1.7882877) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(-0.34298204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(2.2516001) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-0.92372957) q[3];
sx q[3];
rz(-0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(-2.529378) q[0];
rz(1.3587492) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(0.30119687) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33443794) q[0];
sx q[0];
rz(-0.2858735) q[0];
sx q[0];
rz(1.8854333) q[0];
rz(1.1015755) q[2];
sx q[2];
rz(-0.71327268) q[2];
sx q[2];
rz(2.3281946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14872486) q[1];
sx q[1];
rz(-2.5588256) q[1];
sx q[1];
rz(0.73281835) q[1];
rz(-2.4047818) q[3];
sx q[3];
rz(-0.34837803) q[3];
sx q[3];
rz(2.762037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(-0.6855489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-2.8177596) q[0];
rz(-0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-0.30276611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1708831) q[0];
sx q[0];
rz(-1.8237178) q[0];
sx q[0];
rz(-2.80632) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5220171) q[2];
sx q[2];
rz(-1.80772) q[2];
sx q[2];
rz(-2.0120442) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1582022) q[1];
sx q[1];
rz(-1.2081971) q[1];
sx q[1];
rz(2.1561978) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8852194) q[3];
sx q[3];
rz(-1.1629761) q[3];
sx q[3];
rz(0.90424171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6451463) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-0.56078625) q[2];
rz(1.0049817) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(2.0463478) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305369) q[0];
sx q[0];
rz(-2.472214) q[0];
sx q[0];
rz(-0.74991599) q[0];
rz(-2.7610682) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(-0.97533018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98816865) q[0];
sx q[0];
rz(-1.0215217) q[0];
sx q[0];
rz(2.2934521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8504233) q[2];
sx q[2];
rz(-1.5131009) q[2];
sx q[2];
rz(1.1448432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5228032) q[1];
sx q[1];
rz(-2.0674043) q[1];
sx q[1];
rz(1.9267531) q[1];
x q[2];
rz(1.2694025) q[3];
sx q[3];
rz(-1.7835254) q[3];
sx q[3];
rz(-1.0476867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19899496) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(3.0723451) q[3];
sx q[3];
rz(-1.4360177) q[3];
sx q[3];
rz(2.2190905) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0104495) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(-0.3057873) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(-1.1526795) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.688153) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(-2.3896609) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20736097) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(-2.7864252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0704502) q[1];
sx q[1];
rz(-1.8952888) q[1];
sx q[1];
rz(1.9981415) q[1];
x q[2];
rz(1.5074499) q[3];
sx q[3];
rz(-1.5557441) q[3];
sx q[3];
rz(1.2069595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5311188) q[2];
sx q[2];
rz(-2.0788772) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65344812) q[0];
sx q[0];
rz(-1.6076037) q[0];
sx q[0];
rz(-1.6758767) q[0];
rz(2.0987971) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(-2.6063812) q[2];
sx q[2];
rz(-1.3387398) q[2];
sx q[2];
rz(-2.6571318) q[2];
rz(0.10877175) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
