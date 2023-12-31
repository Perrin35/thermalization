OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(4.1783279) q[1];
sx q[1];
rz(8.7611603) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9288869) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(-1.5443718) q[0];
rz(-0.3590091) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(2.5100978) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44890468) q[1];
sx q[1];
rz(-1.9069888) q[1];
sx q[1];
rz(-2.8640139) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1575559) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(1.2260431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(2.7093637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75823443) q[2];
sx q[2];
rz(-0.97823921) q[2];
sx q[2];
rz(-2.638608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.023149816) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(1.5688194) q[1];
x q[2];
rz(-1.5561043) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-2.1121315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7824161) q[0];
sx q[0];
rz(-2.4826907) q[0];
sx q[0];
rz(1.2078148) q[0];
rz(0.6497518) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(2.8284555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.894695) q[1];
sx q[1];
rz(-0.20935911) q[1];
sx q[1];
rz(-0.73114242) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8221178) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(0.6832173) q[0];
rz(1.707294) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(2.0493281) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.013853156) q[1];
sx q[1];
rz(-1.2989559) q[1];
sx q[1];
rz(-2.4747162) q[1];
rz(-pi) q[2];
rz(0.32228542) q[3];
sx q[3];
rz(-2.5282113) q[3];
sx q[3];
rz(1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.8466922) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(0.98658371) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.5198583) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(-2.246726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1631854) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(1.0582256) q[0];
x q[1];
rz(0.61770265) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(-0.4862116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-1.1333229) q[1];
sx q[1];
rz(0.84769627) q[1];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(-2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-0.46498743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156292) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(-2.7727491) q[0];
rz(-0.34157413) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(2.9937033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8919233) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(-3.0628672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6316846) q[0];
sx q[0];
rz(-1.8242867) q[0];
sx q[0];
rz(-1.7476728) q[0];
rz(0.51257001) q[2];
sx q[2];
rz(-1.9480431) q[2];
sx q[2];
rz(0.23114983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5739463) q[1];
sx q[1];
rz(-1.5168377) q[1];
sx q[1];
rz(-1.3957363) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0701253) q[3];
sx q[3];
rz(-1.0896177) q[3];
sx q[3];
rz(-1.037998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.069313958) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-2.8930194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45390689) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.2841671) q[0];
rz(-pi) q[1];
rz(-1.0609264) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(-0.61444297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-0.81070886) q[1];
sx q[1];
rz(-0.87358012) q[1];
x q[2];
rz(-1.400984) q[3];
sx q[3];
rz(-1.6465934) q[3];
sx q[3];
rz(-2.4236987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.0193664) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(3.0526572) q[0];
rz(-0.73762383) q[2];
sx q[2];
rz(-0.80543033) q[2];
sx q[2];
rz(1.9464303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9642155) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(1.7978976) q[1];
rz(2.7301844) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5426853) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.8803966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(-2.7306042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032632685) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(-1.061071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1120158) q[1];
sx q[1];
rz(-2.638991) q[1];
sx q[1];
rz(0.13346787) q[1];
x q[2];
rz(-2.4576549) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(-2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(3.1153395) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(1.431987) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
