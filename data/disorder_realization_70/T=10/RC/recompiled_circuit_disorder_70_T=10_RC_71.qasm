OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4579826) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(-2.3258414) q[0];
rz(-1.6526821) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(1.3002849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12373911) q[1];
sx q[1];
rz(-0.98736073) q[1];
sx q[1];
rz(0.53637335) q[1];
rz(1.5473458) q[3];
sx q[3];
rz(-1.970885) q[3];
sx q[3];
rz(1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76925812) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.9702205) q[0];
rz(-pi) q[1];
rz(-3.0304893) q[2];
sx q[2];
rz(-1.9176033) q[2];
sx q[2];
rz(2.8135516) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-1.1138798) q[1];
sx q[1];
rz(-1.7060075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(0.58829266) q[2];
rz(0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(0.75769889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(-2.589059) q[0];
rz(-pi) q[1];
rz(1.4104841) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(2.2018873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0231087) q[1];
sx q[1];
rz(-2.1493836) q[1];
sx q[1];
rz(-2.5793377) q[1];
rz(0.14962872) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(-1.7865208) q[0];
x q[1];
rz(-0.51283299) q[2];
sx q[2];
rz(-0.83216681) q[2];
sx q[2];
rz(1.2733449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36973876) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(-1.5727444) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10947157) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0230334) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(2.8708007) q[0];
rz(1.1930824) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(1.2010241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5711172) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(-1.2586602) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91315956) q[3];
sx q[3];
rz(-2.1214161) q[3];
sx q[3];
rz(3.1331568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15121962) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(1.3034526) q[0];
rz(0.10411711) q[2];
sx q[2];
rz(-1.8407485) q[2];
sx q[2];
rz(-1.6402668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4420538) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-2.8738408) q[1];
x q[2];
rz(1.7536229) q[3];
sx q[3];
rz(-1.8086686) q[3];
sx q[3];
rz(0.49926234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.2833387) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.3791929) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.798511) q[0];
sx q[0];
rz(-1.3268952) q[0];
sx q[0];
rz(-3.0573465) q[0];
rz(-pi) q[1];
rz(1.9786644) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-2.5025764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52723253) q[1];
sx q[1];
rz(-0.34253201) q[1];
sx q[1];
rz(2.637898) q[1];
x q[2];
rz(-0.93382436) q[3];
sx q[3];
rz(-2.056042) q[3];
sx q[3];
rz(2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.2932628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9613567) q[0];
sx q[0];
rz(-1.1399674) q[0];
sx q[0];
rz(-2.9234617) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5404732) q[2];
sx q[2];
rz(-1.4919315) q[2];
sx q[2];
rz(0.89389801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1845491) q[1];
sx q[1];
rz(-2.7126985) q[1];
sx q[1];
rz(1.9098319) q[1];
x q[2];
rz(0.39956283) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-0.012133908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9171137) q[0];
sx q[0];
rz(-2.0600442) q[0];
sx q[0];
rz(3.1130303) q[0];
rz(-pi) q[1];
x q[1];
rz(1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(0.21208866) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3661256) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(2.8833564) q[1];
rz(-pi) q[2];
rz(-1.6279531) q[3];
sx q[3];
rz(-0.74712979) q[3];
sx q[3];
rz(2.3627918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(2.0690209) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627535) q[0];
sx q[0];
rz(-0.05719962) q[0];
sx q[0];
rz(-2.5662533) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7634002) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(0.80519245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6348833) q[1];
sx q[1];
rz(-2.2193529) q[1];
sx q[1];
rz(-1.3193921) q[1];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(-1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.4762896) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(1.6499741) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(2.9885837) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];