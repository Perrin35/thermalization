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
rz(-1.5795213) q[0];
sx q[0];
rz(6.4223839) q[0];
sx q[0];
rz(9.7425707) q[0];
rz(-2.3843482) q[1];
sx q[1];
rz(4.7321893) q[1];
sx q[1];
rz(7.7319747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4690404) q[0];
sx q[0];
rz(-1.7069478) q[0];
sx q[0];
rz(0.37369136) q[0];
rz(-0.4223675) q[2];
sx q[2];
rz(-2.7701391) q[2];
sx q[2];
rz(-1.708622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9592465) q[1];
sx q[1];
rz(-1.1350147) q[1];
sx q[1];
rz(-2.8027727) q[1];
rz(-1.1640383) q[3];
sx q[3];
rz(-1.6534605) q[3];
sx q[3];
rz(0.57549046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.455787) q[2];
sx q[2];
rz(-1.1263584) q[2];
sx q[2];
rz(-2.8796999) q[2];
rz(-2.453878) q[3];
sx q[3];
rz(-0.7656289) q[3];
sx q[3];
rz(-2.8253637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9460816) q[0];
sx q[0];
rz(-1.2844149) q[0];
sx q[0];
rz(0.84445697) q[0];
rz(-0.63956815) q[1];
sx q[1];
rz(-0.25194672) q[1];
sx q[1];
rz(1.8738939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4019011) q[0];
sx q[0];
rz(-2.2595355) q[0];
sx q[0];
rz(-2.9852377) q[0];
rz(-1.4136281) q[2];
sx q[2];
rz(-0.90051046) q[2];
sx q[2];
rz(-1.8994191) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87593791) q[1];
sx q[1];
rz(-2.3892683) q[1];
sx q[1];
rz(3.1347514) q[1];
rz(-pi) q[2];
rz(0.84898013) q[3];
sx q[3];
rz(-2.0343307) q[3];
sx q[3];
rz(-3.0914538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0212448) q[2];
sx q[2];
rz(-0.52560386) q[2];
sx q[2];
rz(-0.28908238) q[2];
rz(0.31288475) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(0.59278929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.4631571) q[0];
sx q[0];
rz(-2.1364697) q[0];
sx q[0];
rz(2.4682755) q[0];
rz(-2.8939269) q[1];
sx q[1];
rz(-2.1238056) q[1];
sx q[1];
rz(2.8376875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0756098) q[0];
sx q[0];
rz(-2.3572266) q[0];
sx q[0];
rz(0.83215587) q[0];
x q[1];
rz(-0.84776216) q[2];
sx q[2];
rz(-2.1925732) q[2];
sx q[2];
rz(1.1439117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60982882) q[1];
sx q[1];
rz(-1.4312688) q[1];
sx q[1];
rz(-0.05910766) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0237977) q[3];
sx q[3];
rz(-0.46110982) q[3];
sx q[3];
rz(-2.5557704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5690545) q[2];
sx q[2];
rz(-1.3596386) q[2];
sx q[2];
rz(0.016810091) q[2];
rz(-2.861764) q[3];
sx q[3];
rz(-0.40585104) q[3];
sx q[3];
rz(-2.5376937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0957606) q[0];
sx q[0];
rz(-1.1149167) q[0];
sx q[0];
rz(-1.8810077) q[0];
rz(-0.49651217) q[1];
sx q[1];
rz(-1.1505726) q[1];
sx q[1];
rz(1.5042492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53816089) q[0];
sx q[0];
rz(-1.0350845) q[0];
sx q[0];
rz(-2.6441022) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54199647) q[2];
sx q[2];
rz(-1.7769939) q[2];
sx q[2];
rz(-1.3958901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86122266) q[1];
sx q[1];
rz(-0.56764738) q[1];
sx q[1];
rz(1.2595792) q[1];
x q[2];
rz(2.6332985) q[3];
sx q[3];
rz(-2.5807305) q[3];
sx q[3];
rz(1.9702553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.090195) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(-0.16057333) q[2];
rz(1.4422902) q[3];
sx q[3];
rz(-1.5211886) q[3];
sx q[3];
rz(3.1062533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8869239) q[0];
sx q[0];
rz(-1.2441607) q[0];
sx q[0];
rz(0.96920335) q[0];
rz(1.8907549) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(2.9108237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7353775) q[0];
sx q[0];
rz(-0.036140291) q[0];
sx q[0];
rz(-1.4213495) q[0];
x q[1];
rz(-0.0051813263) q[2];
sx q[2];
rz(-1.1924414) q[2];
sx q[2];
rz(-2.3145702) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88420743) q[1];
sx q[1];
rz(-1.1697959) q[1];
sx q[1];
rz(1.0610508) q[1];
rz(-pi) q[2];
rz(2.4668815) q[3];
sx q[3];
rz(-0.58739788) q[3];
sx q[3];
rz(-0.34492465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9427003) q[2];
sx q[2];
rz(-2.5653699) q[2];
sx q[2];
rz(-0.25299859) q[2];
rz(1.3544072) q[3];
sx q[3];
rz(-2.0372882) q[3];
sx q[3];
rz(-1.3033006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.323728) q[0];
sx q[0];
rz(-0.67960328) q[0];
sx q[0];
rz(0.060081765) q[0];
rz(-0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(0.047860535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353226) q[0];
sx q[0];
rz(-1.5859134) q[0];
sx q[0];
rz(1.478805) q[0];
rz(-0.069959241) q[2];
sx q[2];
rz(-0.20393755) q[2];
sx q[2];
rz(-2.4318757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3678494) q[1];
sx q[1];
rz(-1.3700822) q[1];
sx q[1];
rz(1.9108921) q[1];
rz(-pi) q[2];
rz(2.8196857) q[3];
sx q[3];
rz(-1.6471146) q[3];
sx q[3];
rz(-2.8907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(-0.16727373) q[2];
rz(-0.078992756) q[3];
sx q[3];
rz(-1.1827712) q[3];
sx q[3];
rz(-2.3110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1792184) q[0];
sx q[0];
rz(-3.0308864) q[0];
sx q[0];
rz(0.12744823) q[0];
rz(-2.2178862) q[1];
sx q[1];
rz(-2.5741003) q[1];
sx q[1];
rz(-2.4097402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6293242) q[0];
sx q[0];
rz(-1.716181) q[0];
sx q[0];
rz(0.54319197) q[0];
rz(-0.42917128) q[2];
sx q[2];
rz(-1.8597684) q[2];
sx q[2];
rz(1.2850375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7596563) q[1];
sx q[1];
rz(-0.77512533) q[1];
sx q[1];
rz(-1.0785021) q[1];
rz(0.50639373) q[3];
sx q[3];
rz(-0.78575883) q[3];
sx q[3];
rz(0.23791152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38714108) q[2];
sx q[2];
rz(-1.0820505) q[2];
sx q[2];
rz(-0.88073909) q[2];
rz(-0.292101) q[3];
sx q[3];
rz(-2.5160774) q[3];
sx q[3];
rz(-0.97091278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(-2.3604426) q[0];
rz(1.7225522) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(0.17793812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3969361) q[0];
sx q[0];
rz(-1.8783924) q[0];
sx q[0];
rz(1.9768977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4660999) q[2];
sx q[2];
rz(-1.0917821) q[2];
sx q[2];
rz(-1.9164849) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3493578) q[1];
sx q[1];
rz(-0.48949896) q[1];
sx q[1];
rz(1.446318) q[1];
rz(-0.61011647) q[3];
sx q[3];
rz(-0.31829208) q[3];
sx q[3];
rz(1.8221318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6302744) q[2];
sx q[2];
rz(-1.4354458) q[2];
sx q[2];
rz(1.1159631) q[2];
rz(-2.337194) q[3];
sx q[3];
rz(-0.65785995) q[3];
sx q[3];
rz(2.4120954) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6804752) q[0];
sx q[0];
rz(-2.456433) q[0];
sx q[0];
rz(-2.6035736) q[0];
rz(0.33942014) q[1];
sx q[1];
rz(-1.5433886) q[1];
sx q[1];
rz(0.0433878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8755084) q[0];
sx q[0];
rz(-1.5055297) q[0];
sx q[0];
rz(1.4738183) q[0];
rz(-0.15114637) q[2];
sx q[2];
rz(-2.1160963) q[2];
sx q[2];
rz(2.6628833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1307357) q[1];
sx q[1];
rz(-0.82243516) q[1];
sx q[1];
rz(-1.462241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9566986) q[3];
sx q[3];
rz(-1.0066972) q[3];
sx q[3];
rz(-2.7581839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.033120774) q[2];
sx q[2];
rz(-0.35999173) q[2];
sx q[2];
rz(1.6610422) q[2];
rz(2.8451653) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(0.75291434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66210371) q[0];
sx q[0];
rz(-0.39283735) q[0];
sx q[0];
rz(-1.9774849) q[0];
rz(-2.8083943) q[1];
sx q[1];
rz(-1.6644128) q[1];
sx q[1];
rz(0.747917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532497) q[0];
sx q[0];
rz(-2.1624301) q[0];
sx q[0];
rz(1.9337898) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93846069) q[2];
sx q[2];
rz(-1.7424287) q[2];
sx q[2];
rz(2.8454663) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6743636) q[1];
sx q[1];
rz(-1.5780431) q[1];
sx q[1];
rz(-0.054673043) q[1];
rz(-0.87769784) q[3];
sx q[3];
rz(-1.5355645) q[3];
sx q[3];
rz(0.30634634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36772874) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(-1.9749953) q[2];
rz(-1.5424607) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(-2.0994073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3508956) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(-0.028388609) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(-2.1688228) q[2];
sx q[2];
rz(-1.1373873) q[2];
sx q[2];
rz(1.6349229) q[2];
rz(0.52846626) q[3];
sx q[3];
rz(-2.2615643) q[3];
sx q[3];
rz(-0.27558358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
