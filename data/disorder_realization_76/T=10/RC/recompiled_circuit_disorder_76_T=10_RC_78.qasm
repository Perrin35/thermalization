OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(-0.85715357) q[0];
sx q[0];
rz(0.13248086) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24774691) q[0];
sx q[0];
rz(-2.6129299) q[0];
sx q[0];
rz(-1.371944) q[0];
rz(-pi) q[1];
rz(-2.4970826) q[2];
sx q[2];
rz(-1.95382) q[2];
sx q[2];
rz(2.2565947) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8707667) q[1];
sx q[1];
rz(-1.642879) q[1];
sx q[1];
rz(1.4429528) q[1];
rz(-pi) q[2];
rz(0.2487189) q[3];
sx q[3];
rz(-2.0318444) q[3];
sx q[3];
rz(-2.98364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(-1.6202554) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(-0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(1.8923627) q[0];
rz(0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55727977) q[0];
sx q[0];
rz(-1.6788388) q[0];
sx q[0];
rz(2.8312107) q[0];
x q[1];
rz(2.850769) q[2];
sx q[2];
rz(-2.4829156) q[2];
sx q[2];
rz(0.80292279) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3851871) q[1];
sx q[1];
rz(-0.88419534) q[1];
sx q[1];
rz(0.7505721) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35238102) q[3];
sx q[3];
rz(-2.0003194) q[3];
sx q[3];
rz(0.22892117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(-0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648401) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(2.3679249) q[0];
rz(-0.0013008612) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(-0.032827854) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2680227) q[0];
sx q[0];
rz(-0.39186726) q[0];
sx q[0];
rz(-0.64143945) q[0];
rz(-2.122934) q[2];
sx q[2];
rz(-1.3214006) q[2];
sx q[2];
rz(-0.92555911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82959158) q[1];
sx q[1];
rz(-2.3288576) q[1];
sx q[1];
rz(0.22711046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0333943) q[3];
sx q[3];
rz(-0.61363797) q[3];
sx q[3];
rz(-2.9463241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(-0.2335877) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(-0.77082005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642826) q[0];
sx q[0];
rz(-1.5872123) q[0];
sx q[0];
rz(1.5903227) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0609444) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(0.25445081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0499038) q[1];
sx q[1];
rz(-0.65444512) q[1];
sx q[1];
rz(0.62033886) q[1];
rz(-pi) q[2];
rz(1.3464438) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(-2.2546774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(0.7129933) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(-0.094093181) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(0.17000155) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029179137) q[0];
sx q[0];
rz(-2.3110483) q[0];
sx q[0];
rz(-1.889617) q[0];
rz(1.1414358) q[2];
sx q[2];
rz(-2.6300207) q[2];
sx q[2];
rz(2.9432952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4449094) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(-2.6615104) q[1];
rz(-pi) q[2];
rz(-0.65808987) q[3];
sx q[3];
rz(-2.1562025) q[3];
sx q[3];
rz(-2.493849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(1.4040995) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(0.45853841) q[0];
rz(0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(0.68516723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72502575) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(-0.91992232) q[0];
x q[1];
rz(-0.14641996) q[2];
sx q[2];
rz(-2.358987) q[2];
sx q[2];
rz(-2.866982) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7428776) q[1];
sx q[1];
rz(-1.0814953) q[1];
sx q[1];
rz(2.2687885) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1232883) q[3];
sx q[3];
rz(-0.98494512) q[3];
sx q[3];
rz(2.626112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3926065) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(-0.72934735) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.9940631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097737) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(-0.41418196) q[0];
rz(-pi) q[1];
rz(-2.0562901) q[2];
sx q[2];
rz(-0.69316961) q[2];
sx q[2];
rz(-1.327141) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3299678) q[1];
sx q[1];
rz(-2.6553272) q[1];
sx q[1];
rz(0.087841308) q[1];
rz(-pi) q[2];
rz(2.6461584) q[3];
sx q[3];
rz(-1.9154275) q[3];
sx q[3];
rz(-0.37638327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-0.914004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4115899) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(1.0986064) q[0];
rz(-pi) q[1];
x q[1];
rz(2.987791) q[2];
sx q[2];
rz(-0.6066583) q[2];
sx q[2];
rz(-0.80468824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1917845) q[1];
sx q[1];
rz(-2.5439918) q[1];
sx q[1];
rz(1.530184) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(-0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.297696) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(1.3893611) q[0];
rz(1.5147491) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(1.0983889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669945) q[0];
sx q[0];
rz(-0.54531389) q[0];
sx q[0];
rz(1.1101515) q[0];
rz(-0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(2.2965477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92022773) q[1];
sx q[1];
rz(-0.5792633) q[1];
sx q[1];
rz(-3.1406162) q[1];
x q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(-0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0825723) q[0];
sx q[0];
rz(-1.1963545) q[0];
sx q[0];
rz(-3.0890205) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8877108) q[2];
sx q[2];
rz(-0.66714087) q[2];
sx q[2];
rz(2.4184879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(0.44209977) q[1];
x q[2];
rz(-2.4217989) q[3];
sx q[3];
rz(-1.8344973) q[3];
sx q[3];
rz(-1.3849474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24370596) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(-1.0591327) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(0.89390755) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(-0.25390608) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(0.72348307) q[2];
sx q[2];
rz(-1.537848) q[2];
sx q[2];
rz(1.6707735) q[2];
rz(-0.068594882) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
