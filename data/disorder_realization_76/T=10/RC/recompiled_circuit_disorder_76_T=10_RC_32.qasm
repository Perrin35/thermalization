OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(3.0091118) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(1.7696487) q[0];
rz(2.5508467) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(-0.22434805) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8707667) q[1];
sx q[1];
rz(-1.4987136) q[1];
sx q[1];
rz(1.4429528) q[1];
rz(2.8928738) q[3];
sx q[3];
rz(-1.1097483) q[3];
sx q[3];
rz(0.15795262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0341558) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(1.6050603) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261616) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97894325) q[0];
sx q[0];
rz(-1.2622841) q[0];
sx q[0];
rz(-1.4573775) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3524019) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(-0.44109694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3483352) q[1];
sx q[1];
rz(-2.1264592) q[1];
sx q[1];
rz(2.4130915) q[1];
rz(-pi) q[2];
rz(2.7892116) q[3];
sx q[3];
rz(-2.0003194) q[3];
sx q[3];
rz(2.9126715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4124174) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(2.5382606) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(0.032827854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0484682) q[0];
sx q[0];
rz(-1.801352) q[0];
sx q[0];
rz(-2.8218517) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0186586) q[2];
sx q[2];
rz(-1.3214006) q[2];
sx q[2];
rz(2.2160335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(-2.3418531) q[1];
rz(-pi) q[2];
rz(0.34605108) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924776) q[0];
sx q[0];
rz(-0.025509398) q[0];
sx q[0];
rz(-0.87164469) q[0];
rz(1.3768251) q[2];
sx q[2];
rz(-2.6435404) q[2];
sx q[2];
rz(1.9961403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3168287) q[1];
sx q[1];
rz(-1.0526122) q[1];
sx q[1];
rz(1.1512685) q[1];
rz(-pi) q[2];
rz(-1.3464438) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(2.2546774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9757441) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029179137) q[0];
sx q[0];
rz(-2.3110483) q[0];
sx q[0];
rz(-1.889617) q[0];
x q[1];
rz(0.22959392) q[2];
sx q[2];
rz(-2.0320971) q[2];
sx q[2];
rz(-0.68179828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4449094) q[1];
sx q[1];
rz(-2.141436) q[1];
sx q[1];
rz(-2.6615104) q[1];
x q[2];
rz(-2.2682297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(-1.8148592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(1.5480301) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858522) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(2.4564254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348649) q[0];
sx q[0];
rz(-0.92175882) q[0];
sx q[0];
rz(0.087260212) q[0];
rz(-pi) q[1];
rz(0.77722705) q[2];
sx q[2];
rz(-1.4677375) q[2];
sx q[2];
rz(1.7412141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8243858) q[1];
sx q[1];
rz(-2.3134391) q[1];
sx q[1];
rz(-2.2627027) q[1];
rz(-pi) q[2];
rz(-0.98486793) q[3];
sx q[3];
rz(-1.586048) q[3];
sx q[3];
rz(-2.0963984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012506164) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(-2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.9940631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80619752) q[0];
sx q[0];
rz(-1.9837539) q[0];
sx q[0];
rz(-1.6523244) q[0];
rz(-2.0562901) q[2];
sx q[2];
rz(-0.69316961) q[2];
sx q[2];
rz(1.8144516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9784669) q[1];
sx q[1];
rz(-1.5297869) q[1];
sx q[1];
rz(-0.48467111) q[1];
x q[2];
rz(1.958193) q[3];
sx q[3];
rz(-2.0347188) q[3];
sx q[3];
rz(-2.1277609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(-2.0429862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15380165) q[2];
sx q[2];
rz(-0.6066583) q[2];
sx q[2];
rz(-0.80468824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1426705) q[1];
sx q[1];
rz(-2.167836) q[1];
sx q[1];
rz(0.02762694) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0141482) q[3];
sx q[3];
rz(-2.017052) q[3];
sx q[3];
rz(-0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(-1.297696) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(-1.5147491) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(1.0983889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3745981) q[0];
sx q[0];
rz(-0.54531389) q[0];
sx q[0];
rz(-1.1101515) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(-2.2965477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64975148) q[1];
sx q[1];
rz(-1.5713308) q[1];
sx q[1];
rz(2.5623296) q[1];
rz(-2.5536355) q[3];
sx q[3];
rz(-0.30246099) q[3];
sx q[3];
rz(1.0126225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(-1.3174723) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53101978) q[0];
sx q[0];
rz(-1.5218698) q[0];
sx q[0];
rz(-1.945709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8877108) q[2];
sx q[2];
rz(-0.66714087) q[2];
sx q[2];
rz(0.72310477) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(2.6994929) q[1];
x q[2];
rz(0.38871308) q[3];
sx q[3];
rz(-0.75838381) q[3];
sx q[3];
rz(0.10314108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(-2.8876866) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(0.72348307) q[2];
sx q[2];
rz(-1.537848) q[2];
sx q[2];
rz(1.6707735) q[2];
rz(-3.0729978) q[3];
sx q[3];
rz(-0.98496901) q[3];
sx q[3];
rz(1.2377644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
