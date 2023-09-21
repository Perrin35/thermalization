OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(1.5443718) q[0];
x q[1];
rz(0.3590091) q[2];
sx q[2];
rz(-0.81072545) q[2];
sx q[2];
rz(2.5100978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8780958) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(-0.90579512) q[1];
rz(-pi) q[2];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7694089) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(2.7093637) q[0];
x q[1];
rz(0.77483564) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(-1.6006084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1204685) q[1];
sx q[1];
rz(-1.7906584) q[1];
sx q[1];
rz(3.1411509) q[1];
rz(-pi) q[2];
rz(2.1434104) q[3];
sx q[3];
rz(-1.5628353) q[3];
sx q[3];
rz(1.5496467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(-0.94421454) q[0];
x q[1];
rz(2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(-1.5497269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50476915) q[1];
sx q[1];
rz(-1.7261337) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.401628) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(1.9268074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8221178) q[0];
sx q[0];
rz(-1.7135156) q[0];
sx q[0];
rz(2.4583754) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0558673) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(-0.9312219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.792946) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(-1.2299041) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5528615) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(-0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.8466922) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-0.98658371) q[3];
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
rz(-2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2762404) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-2.4549237) q[0];
rz(2.4159555) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(2.5051136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4836854) q[1];
sx q[1];
rz(-2.2135418) q[1];
sx q[1];
rz(-0.55773463) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1154247) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(0.22931306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.616509) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(2.6766052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156292) q[0];
sx q[0];
rz(-2.0032126) q[0];
sx q[0];
rz(-0.36884357) q[0];
rz(-pi) q[1];
rz(0.34157413) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(-2.9937033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28593674) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(-3.082285) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(3.0628672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(2.6965551) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10931817) q[0];
sx q[0];
rz(-2.8335857) q[0];
sx q[0];
rz(-0.59662915) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4620352) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(-1.9192413) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.012688795) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(-3.0867982) q[1];
x q[2];
rz(-2.0701253) q[3];
sx q[3];
rz(-2.0519749) q[3];
sx q[3];
rz(-1.037998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(1.0726661) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45390689) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.8574255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82069223) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0174745) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(2.2494621) q[1];
x q[2];
rz(-3.0646938) q[3];
sx q[3];
rz(-1.7401164) q[3];
sx q[3];
rz(0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.3185906) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5134207) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(1.1286939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4852072) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(-0.93751794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88398891) q[1];
sx q[1];
rz(-2.249243) q[1];
sx q[1];
rz(0.18737327) q[1];
x q[2];
rz(0.41140822) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(-2.0421162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7302154) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.8803966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48860088) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(-2.7306042) q[0];
x q[1];
rz(2.5433259) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(2.604904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0295769) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(-3.0081248) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4576549) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(-0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(1.6115887) q[2];
sx q[2];
rz(-0.57208021) q[2];
sx q[2];
rz(-0.013442599) q[2];
rz(-2.5504997) q[3];
sx q[3];
rz(-1.686284) q[3];
sx q[3];
rz(-1.4808663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
