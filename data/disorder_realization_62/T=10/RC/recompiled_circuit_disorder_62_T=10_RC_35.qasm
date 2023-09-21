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
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463319) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(-0.46121009) q[0];
rz(1.2167591) q[2];
sx q[2];
rz(-0.82497043) q[2];
sx q[2];
rz(2.0113457) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44890468) q[1];
sx q[1];
rz(-1.9069888) q[1];
sx q[1];
rz(-0.2775788) q[1];
rz(-pi) q[2];
rz(3.033028) q[3];
sx q[3];
rz(-1.8130455) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(0.051068548) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(1.9155496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2873043) q[0];
sx q[0];
rz(-1.1278296) q[0];
sx q[0];
rz(-1.3501549) q[0];
rz(0.82310279) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(2.5603103) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5495758) q[1];
sx q[1];
rz(-1.5712275) q[1];
sx q[1];
rz(1.7906584) q[1];
x q[2];
rz(1.5854884) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80699608) q[0];
sx q[0];
rz(-0.96141978) q[0];
sx q[0];
rz(-0.2683123) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2842032) q[2];
sx q[2];
rz(-2.2006052) q[2];
sx q[2];
rz(2.0558002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0440812) q[1];
sx q[1];
rz(-1.4315726) q[1];
sx q[1];
rz(0.15686762) q[1];
x q[2];
rz(2.401628) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(-1.9268074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.625659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-0.69561361) q[0];
sx q[0];
rz(0.22380933) q[0];
rz(-1.4342986) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(2.0493281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.792946) q[1];
sx q[1];
rz(-2.209084) q[1];
sx q[1];
rz(-1.9116886) q[1];
x q[2];
rz(-0.32228542) q[3];
sx q[3];
rz(-2.5282113) q[3];
sx q[3];
rz(-1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(-1.5198583) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1631854) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(2.083367) q[0];
rz(-0.61770265) q[2];
sx q[2];
rz(-2.3139944) q[2];
sx q[2];
rz(-0.4862116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.289031) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(2.186071) q[1];
rz(3.1154247) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(0.46498743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573479) q[0];
sx q[0];
rz(-1.9042958) q[0];
sx q[0];
rz(-2.0302982) q[0];
x q[1];
rz(0.20435135) q[2];
sx q[2];
rz(-0.34826476) q[2];
sx q[2];
rz(-1.5262926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28593674) q[1];
sx q[1];
rz(-2.2054513) q[1];
sx q[1];
rz(0.059307701) q[1];
x q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(3.0628672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359081) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(-2.8842584) q[0];
rz(-0.67955741) q[2];
sx q[2];
rz(-0.62629269) q[2];
sx q[2];
rz(-1.9192413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1289039) q[1];
sx q[1];
rz(-1.745599) q[1];
sx q[1];
rz(0.054794475) q[1];
x q[2];
rz(0.74219269) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(-1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.5163039) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
sx q[2];
rz(pi/2) q[2];
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
rz(2.3209004) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2182902) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(0.87358012) q[1];
rz(1.400984) q[3];
sx q[3];
rz(-1.6465934) q[3];
sx q[3];
rz(2.4236987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2934072) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(-1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(-2.0128987) q[0];
x q[1];
rz(-0.65638541) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(0.93751794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56837592) q[1];
sx q[1];
rz(-1.7163367) q[1];
sx q[1];
rz(2.2578866) q[1];
rz(1.9570458) q[3];
sx q[3];
rz(-2.2830314) q[3];
sx q[3];
rz(-1.6561179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3189145) q[0];
sx q[0];
rz(-0.42352522) q[0];
sx q[0];
rz(-2.8828997) q[0];
rz(-0.59826675) q[2];
sx q[2];
rz(-1.5524128) q[2];
sx q[2];
rz(0.53668864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.48325) q[1];
sx q[1];
rz(-1.634942) q[1];
sx q[1];
rz(-2.642753) q[1];
rz(-2.4576549) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(-0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
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
rz(0.026253168) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(0.20523397) q[3];
sx q[3];
rz(-2.5406465) q[3];
sx q[3];
rz(3.061486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];