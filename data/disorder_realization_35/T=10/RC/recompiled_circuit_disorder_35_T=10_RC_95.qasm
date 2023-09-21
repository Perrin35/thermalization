OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5380733) q[0];
sx q[0];
rz(-2.3581714) q[0];
sx q[0];
rz(0.51622434) q[0];
rz(-pi) q[1];
rz(-2.3697882) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(-0.31847218) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.484326) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(-2.3372997) q[1];
rz(-pi) q[2];
rz(2.911525) q[3];
sx q[3];
rz(-0.32031968) q[3];
sx q[3];
rz(2.7473118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(-0.47505501) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(2.1038726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25181928) q[0];
sx q[0];
rz(-0.60476859) q[0];
sx q[0];
rz(0.39321995) q[0];
x q[1];
rz(0.53595397) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(0.11975372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71643752) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(-1.1643216) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.188835) q[3];
sx q[3];
rz(-2.2817094) q[3];
sx q[3];
rz(-3.0666805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(2.5684165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37704913) q[0];
sx q[0];
rz(-1.1408313) q[0];
sx q[0];
rz(3.0717875) q[0];
rz(2.8163221) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(-0.48649597) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0108311) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-0.39949135) q[1];
rz(-pi) q[2];
rz(-1.1099986) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4578611) q[0];
sx q[0];
rz(-2.3944003) q[0];
sx q[0];
rz(-1.9283717) q[0];
x q[1];
rz(-2.0998459) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(-1.0066102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73211654) q[1];
sx q[1];
rz(-0.64142694) q[1];
sx q[1];
rz(-0.17318053) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(-2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-0.0030227946) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.682122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057102324) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(-1.811972) q[0];
x q[1];
rz(2.2541788) q[2];
sx q[2];
rz(-1.9878584) q[2];
sx q[2];
rz(-0.29287698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0046878) q[1];
sx q[1];
rz(-1.3740174) q[1];
sx q[1];
rz(2.911527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8206283) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(0.13866436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16620557) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-0.020394527) q[0];
rz(-pi) q[1];
rz(-1.6476829) q[2];
sx q[2];
rz(-1.7711519) q[2];
sx q[2];
rz(-2.4424057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8530635) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(2.2437375) q[1];
rz(1.9964553) q[3];
sx q[3];
rz(-0.87493757) q[3];
sx q[3];
rz(-1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(0.62430635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(1.2929582) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1820071) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(2.6401273) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.034060409) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(2.2570616) q[1];
rz(-0.94654406) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(3.1012227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-3.0126742) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3710204) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(-1.8380941) q[0];
x q[1];
rz(2.244197) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(-0.71133864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32313777) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(2.6998991) q[1];
rz(-pi) q[2];
rz(1.586732) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(2.074266) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006376) q[0];
sx q[0];
rz(-1.3906286) q[0];
sx q[0];
rz(0.10794497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9753014) q[2];
sx q[2];
rz(-2.1098237) q[2];
sx q[2];
rz(-1.5664958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0061134) q[1];
sx q[1];
rz(-2.0524426) q[1];
sx q[1];
rz(-2.3856132) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25165598) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(-2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71781681) q[0];
sx q[0];
rz(-2.2795838) q[0];
sx q[0];
rz(0.29341673) q[0];
rz(-0.28995138) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(1.314756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5200561) q[1];
sx q[1];
rz(-1.6076419) q[1];
sx q[1];
rz(0.28921339) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5952026) q[3];
sx q[3];
rz(-0.26248172) q[3];
sx q[3];
rz(2.4266092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(-0.98106445) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];