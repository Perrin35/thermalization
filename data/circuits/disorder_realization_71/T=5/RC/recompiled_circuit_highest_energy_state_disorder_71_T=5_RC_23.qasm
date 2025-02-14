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
rz(-2.0105536) q[0];
sx q[0];
rz(-2.5656963) q[0];
sx q[0];
rz(2.0146712) q[0];
rz(2.1908886) q[1];
sx q[1];
rz(3.7428441) q[1];
sx q[1];
rz(9.7584702) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5883049) q[0];
sx q[0];
rz(-1.9899564) q[0];
sx q[0];
rz(2.9774211) q[0];
x q[1];
rz(1.1154217) q[2];
sx q[2];
rz(-1.3297594) q[2];
sx q[2];
rz(-0.068403989) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0036949) q[1];
sx q[1];
rz(-2.7337791) q[1];
sx q[1];
rz(0.56038709) q[1];
rz(1.8435616) q[3];
sx q[3];
rz(-2.354151) q[3];
sx q[3];
rz(-1.8026601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8814964) q[2];
sx q[2];
rz(-1.0755971) q[2];
sx q[2];
rz(2.0913701) q[2];
rz(2.8460734) q[3];
sx q[3];
rz(-0.76885709) q[3];
sx q[3];
rz(1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1693901) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(2.4743359) q[0];
rz(0.15070209) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(2.8299455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.920645) q[0];
sx q[0];
rz(-1.049982) q[0];
sx q[0];
rz(-3.060311) q[0];
rz(-1.1318593) q[2];
sx q[2];
rz(-1.2905741) q[2];
sx q[2];
rz(-2.9228022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91193059) q[1];
sx q[1];
rz(-1.4313102) q[1];
sx q[1];
rz(0.18499495) q[1];
x q[2];
rz(2.7439609) q[3];
sx q[3];
rz(-1.5036229) q[3];
sx q[3];
rz(-0.45715082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47505891) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(2.3739572) q[2];
rz(-2.660699) q[3];
sx q[3];
rz(-2.3700263) q[3];
sx q[3];
rz(-1.5868928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84592205) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(0.77958244) q[0];
rz(0.37173158) q[1];
sx q[1];
rz(-1.0075684) q[1];
sx q[1];
rz(-0.76098162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62255732) q[0];
sx q[0];
rz(-1.8599959) q[0];
sx q[0];
rz(-1.2361121) q[0];
rz(-pi) q[1];
rz(2.089431) q[2];
sx q[2];
rz(-2.385879) q[2];
sx q[2];
rz(1.99988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5058712) q[1];
sx q[1];
rz(-2.2835915) q[1];
sx q[1];
rz(0.21943032) q[1];
rz(-pi) q[2];
rz(-0.0049505135) q[3];
sx q[3];
rz(-0.85270665) q[3];
sx q[3];
rz(2.6867721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80428213) q[2];
sx q[2];
rz(-2.2411942) q[2];
sx q[2];
rz(0.4551746) q[2];
rz(0.69194397) q[3];
sx q[3];
rz(-0.71440905) q[3];
sx q[3];
rz(-0.80879319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96875018) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(-0.2562879) q[0];
rz(-1.434727) q[1];
sx q[1];
rz(-1.7211569) q[1];
sx q[1];
rz(-3.0228379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.502536) q[0];
sx q[0];
rz(-2.0169741) q[0];
sx q[0];
rz(1.3986375) q[0];
rz(-pi) q[1];
rz(1.9740749) q[2];
sx q[2];
rz(-0.31709798) q[2];
sx q[2];
rz(-0.69363475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.465047) q[1];
sx q[1];
rz(-0.80261723) q[1];
sx q[1];
rz(1.8957183) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7440192) q[3];
sx q[3];
rz(-2.3459917) q[3];
sx q[3];
rz(-0.52594409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8541096) q[2];
sx q[2];
rz(-0.27366769) q[2];
sx q[2];
rz(-1.0556833) q[2];
rz(1.4500729) q[3];
sx q[3];
rz(-1.6943211) q[3];
sx q[3];
rz(2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56181041) q[0];
sx q[0];
rz(-0.17794839) q[0];
sx q[0];
rz(-1.6135038) q[0];
rz(-1.9644507) q[1];
sx q[1];
rz(-1.2700932) q[1];
sx q[1];
rz(1.5632163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684473) q[0];
sx q[0];
rz(-0.6296491) q[0];
sx q[0];
rz(-0.11778732) q[0];
x q[1];
rz(-1.8991532) q[2];
sx q[2];
rz(-1.17982) q[2];
sx q[2];
rz(2.7308488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0381752) q[1];
sx q[1];
rz(-0.6486054) q[1];
sx q[1];
rz(-2.9056249) q[1];
rz(-pi) q[2];
rz(-1.5734736) q[3];
sx q[3];
rz(-1.9802046) q[3];
sx q[3];
rz(-2.1684627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13581181) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(1.8394252) q[2];
rz(-0.33995134) q[3];
sx q[3];
rz(-0.8067185) q[3];
sx q[3];
rz(2.4792041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.133701) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.1997461) q[0];
rz(1.7802736) q[1];
sx q[1];
rz(-2.168455) q[1];
sx q[1];
rz(0.57788411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21450522) q[0];
sx q[0];
rz(-1.5643692) q[0];
sx q[0];
rz(-1.584573) q[0];
rz(-1.4125713) q[2];
sx q[2];
rz(-1.998868) q[2];
sx q[2];
rz(-2.9077663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2924345) q[1];
sx q[1];
rz(-0.47579381) q[1];
sx q[1];
rz(-2.8413296) q[1];
rz(1.2718926) q[3];
sx q[3];
rz(-1.546197) q[3];
sx q[3];
rz(-0.92211039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5100688) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(-0.4734545) q[2];
rz(1.2724642) q[3];
sx q[3];
rz(-1.3926287) q[3];
sx q[3];
rz(2.9191391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1943787) q[0];
sx q[0];
rz(-0.65776238) q[0];
sx q[0];
rz(-2.8084602) q[0];
rz(-0.22854742) q[1];
sx q[1];
rz(-1.3401778) q[1];
sx q[1];
rz(-2.5659335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822214) q[0];
sx q[0];
rz(-0.8176119) q[0];
sx q[0];
rz(-2.9961917) q[0];
x q[1];
rz(-2.5533122) q[2];
sx q[2];
rz(-0.55598393) q[2];
sx q[2];
rz(0.54186547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0488494) q[1];
sx q[1];
rz(-1.7705838) q[1];
sx q[1];
rz(-1.9034991) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3153399) q[3];
sx q[3];
rz(-1.1029579) q[3];
sx q[3];
rz(1.6088736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8754742) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(-1.4761338) q[2];
rz(2.0406145) q[3];
sx q[3];
rz(-1.063238) q[3];
sx q[3];
rz(-0.5180009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440893) q[0];
sx q[0];
rz(-0.88633716) q[0];
sx q[0];
rz(-2.8344717) q[0];
rz(0.33044526) q[1];
sx q[1];
rz(-1.5978866) q[1];
sx q[1];
rz(-1.365136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70676409) q[0];
sx q[0];
rz(-1.4894823) q[0];
sx q[0];
rz(1.5887194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34806378) q[2];
sx q[2];
rz(-1.3937147) q[2];
sx q[2];
rz(-0.40906104) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.04330895) q[1];
sx q[1];
rz(-1.829147) q[1];
sx q[1];
rz(-2.0140548) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7264113) q[3];
sx q[3];
rz(-2.3978516) q[3];
sx q[3];
rz(1.9187601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27552989) q[2];
sx q[2];
rz(-2.5538462) q[2];
sx q[2];
rz(-2.8928939) q[2];
rz(1.2134877) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(0.061633751) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520749) q[0];
sx q[0];
rz(-2.2353421) q[0];
sx q[0];
rz(0.686598) q[0];
rz(-1.1129414) q[1];
sx q[1];
rz(-0.80250347) q[1];
sx q[1];
rz(-0.31563219) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7677143) q[0];
sx q[0];
rz(-1.6168103) q[0];
sx q[0];
rz(-1.893928) q[0];
rz(0.2608582) q[2];
sx q[2];
rz(-2.3645698) q[2];
sx q[2];
rz(2.379247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3179847) q[1];
sx q[1];
rz(-2.4725742) q[1];
sx q[1];
rz(1.7229592) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0635438) q[3];
sx q[3];
rz(-2.4025318) q[3];
sx q[3];
rz(-0.8821677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.687872) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(0.68457121) q[2];
rz(1.941393) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(-2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612563) q[0];
sx q[0];
rz(-2.6604524) q[0];
sx q[0];
rz(-0.0048333724) q[0];
rz(-1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(2.8878816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8953573) q[0];
sx q[0];
rz(-2.9140436) q[0];
sx q[0];
rz(1.958311) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6275835) q[2];
sx q[2];
rz(-0.80482641) q[2];
sx q[2];
rz(0.45586205) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.743266) q[1];
sx q[1];
rz(-0.70320819) q[1];
sx q[1];
rz(1.2656256) q[1];
rz(2.8242094) q[3];
sx q[3];
rz(-0.57358984) q[3];
sx q[3];
rz(1.7490179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1634875) q[2];
sx q[2];
rz(-1.9065964) q[2];
sx q[2];
rz(-1.9592436) q[2];
rz(-0.67665082) q[3];
sx q[3];
rz(-0.38656056) q[3];
sx q[3];
rz(1.0138018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61426281) q[0];
sx q[0];
rz(-0.78582055) q[0];
sx q[0];
rz(0.2926122) q[0];
rz(1.0489427) q[1];
sx q[1];
rz(-1.4194149) q[1];
sx q[1];
rz(0.70379757) q[1];
rz(-2.537773) q[2];
sx q[2];
rz(-0.55470822) q[2];
sx q[2];
rz(1.8798238) q[2];
rz(1.3181237) q[3];
sx q[3];
rz(-0.37794401) q[3];
sx q[3];
rz(1.5156802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
