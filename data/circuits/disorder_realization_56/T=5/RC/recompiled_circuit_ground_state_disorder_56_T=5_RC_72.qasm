OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(0.78328744) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(-2.3496871) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053122205) q[0];
sx q[0];
rz(-1.7566009) q[0];
sx q[0];
rz(2.8775999) q[0];
rz(1.3375086) q[2];
sx q[2];
rz(-2.983494) q[2];
sx q[2];
rz(2.5776517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72630771) q[1];
sx q[1];
rz(-2.1139164) q[1];
sx q[1];
rz(-1.607479) q[1];
rz(-pi) q[2];
rz(-0.79444076) q[3];
sx q[3];
rz(-1.1367961) q[3];
sx q[3];
rz(1.7931012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59114328) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(0.24672306) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-0.39768404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(0.9285399) q[0];
rz(-1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(-0.64847747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47005338) q[0];
sx q[0];
rz(-0.82255615) q[0];
sx q[0];
rz(0.7732735) q[0];
x q[1];
rz(-1.8869867) q[2];
sx q[2];
rz(-2.11497) q[2];
sx q[2];
rz(0.7989102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3120756) q[1];
sx q[1];
rz(-1.5863998) q[1];
sx q[1];
rz(0.14194686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8594555) q[3];
sx q[3];
rz(-1.6842173) q[3];
sx q[3];
rz(-2.3577549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4497946) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(-2.2890384) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5071252) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(-2.0563828) q[0];
rz(-0.94404864) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(-1.5012213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0414934) q[0];
sx q[0];
rz(-1.7679119) q[0];
sx q[0];
rz(1.5310578) q[0];
x q[1];
rz(-2.9502421) q[2];
sx q[2];
rz(-1.3406959) q[2];
sx q[2];
rz(-1.5887345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.855401) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(-1.3894677) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7872058) q[3];
sx q[3];
rz(-1.7248271) q[3];
sx q[3];
rz(-2.0311525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40470716) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(-0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(-1.1165718) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214355) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(-1.1965363) q[0];
rz(1.6664956) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(1.9570785) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329342) q[0];
sx q[0];
rz(-2.2581165) q[0];
sx q[0];
rz(1.2615027) q[0];
x q[1];
rz(-2.938561) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(-2.610746) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8912011) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(-2.0865284) q[1];
rz(-pi) q[2];
rz(1.3877901) q[3];
sx q[3];
rz(-1.840045) q[3];
sx q[3];
rz(1.315801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2148332) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(2.7091806) q[3];
sx q[3];
rz(-2.0093446) q[3];
sx q[3];
rz(0.47237083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.072902) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(-0.75621653) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(-2.3888033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469893) q[0];
sx q[0];
rz(-1.5651817) q[0];
sx q[0];
rz(3.1412197) q[0];
rz(-pi) q[1];
rz(3.0514206) q[2];
sx q[2];
rz(-1.3036696) q[2];
sx q[2];
rz(1.4324485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.237287) q[1];
sx q[1];
rz(-0.5798961) q[1];
sx q[1];
rz(-0.19176264) q[1];
x q[2];
rz(-0.54075586) q[3];
sx q[3];
rz(-1.9633368) q[3];
sx q[3];
rz(-0.3730216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14043643) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(0.39503869) q[2];
rz(2.6664074) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(0.37676677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(2.6267701) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(2.8033676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048934919) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(0.2661163) q[0];
x q[1];
rz(-0.80775308) q[2];
sx q[2];
rz(-2.8050426) q[2];
sx q[2];
rz(1.6405662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5023313) q[1];
sx q[1];
rz(-0.0031454589) q[1];
sx q[1];
rz(1.5793213) q[1];
rz(3.0727524) q[3];
sx q[3];
rz(-1.9919812) q[3];
sx q[3];
rz(2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6414791) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8165269) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.3651747) q[1];
sx q[1];
rz(-0.65779984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938134) q[0];
sx q[0];
rz(-1.9792611) q[0];
sx q[0];
rz(2.4373217) q[0];
rz(1.5267685) q[2];
sx q[2];
rz(-2.2775473) q[2];
sx q[2];
rz(2.6296774) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15193394) q[1];
sx q[1];
rz(-0.2865782) q[1];
sx q[1];
rz(-1.6295939) q[1];
x q[2];
rz(-0.21605394) q[3];
sx q[3];
rz(-1.1974058) q[3];
sx q[3];
rz(-1.4114398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49360069) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(-1.4823401) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7928612) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(-0.29801512) q[0];
rz(-1.0385849) q[1];
sx q[1];
rz(-0.54439259) q[1];
sx q[1];
rz(-2.9878152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6385495) q[0];
sx q[0];
rz(-2.8766603) q[0];
sx q[0];
rz(-2.8517836) q[0];
x q[1];
rz(0.80259097) q[2];
sx q[2];
rz(-2.0611412) q[2];
sx q[2];
rz(-2.5732694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41339918) q[1];
sx q[1];
rz(-1.4811885) q[1];
sx q[1];
rz(0.34743584) q[1];
rz(0.12184398) q[3];
sx q[3];
rz(-1.9949556) q[3];
sx q[3];
rz(-2.3736853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0457354) q[2];
sx q[2];
rz(-2.6726674) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(-0.73330283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.34340149) q[0];
sx q[0];
rz(-2.5372086) q[0];
sx q[0];
rz(1.6424302) q[0];
rz(0.54939735) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3651705) q[0];
sx q[0];
rz(-0.93342268) q[0];
sx q[0];
rz(2.754162) q[0];
rz(-pi) q[1];
rz(2.8806503) q[2];
sx q[2];
rz(-1.5847209) q[2];
sx q[2];
rz(1.5557289) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10560606) q[1];
sx q[1];
rz(-2.3422554) q[1];
sx q[1];
rz(1.6677594) q[1];
rz(-2.9467877) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(2.3969802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(-0.086183444) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(-2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8330399) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(0.53681701) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(-0.73582617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9402855) q[0];
sx q[0];
rz(-0.92865151) q[0];
sx q[0];
rz(1.3298195) q[0];
x q[1];
rz(0.64401099) q[2];
sx q[2];
rz(-2.0154872) q[2];
sx q[2];
rz(1.4294238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5246626) q[1];
sx q[1];
rz(-1.7866275) q[1];
sx q[1];
rz(2.4991577) q[1];
rz(-pi) q[2];
rz(2.5352468) q[3];
sx q[3];
rz(-0.56713533) q[3];
sx q[3];
rz(0.66886574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95911038) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(2.518892) q[2];
rz(-0.73838082) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(2.8625989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57711346) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(-0.64361698) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(-1.2014482) q[2];
sx q[2];
rz(-1.7864173) q[2];
sx q[2];
rz(-1.3997072) q[2];
rz(-2.6162061) q[3];
sx q[3];
rz(-2.864884) q[3];
sx q[3];
rz(-1.3340193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
