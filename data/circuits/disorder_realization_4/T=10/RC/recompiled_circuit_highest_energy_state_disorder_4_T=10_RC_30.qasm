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
rz(-1.1308489) q[0];
sx q[0];
rz(-1.3698438) q[0];
sx q[0];
rz(-0.1432336) q[0];
rz(3.1205966) q[1];
sx q[1];
rz(3.723998) q[1];
sx q[1];
rz(5.8082685) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62014318) q[0];
sx q[0];
rz(-3.0472429) q[0];
sx q[0];
rz(-0.63395377) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5449499) q[2];
sx q[2];
rz(-0.95603564) q[2];
sx q[2];
rz(1.3790707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7066267) q[1];
sx q[1];
rz(-2.3449223) q[1];
sx q[1];
rz(-0.15310751) q[1];
rz(-pi) q[2];
rz(-2.0359614) q[3];
sx q[3];
rz(-2.0099927) q[3];
sx q[3];
rz(0.10094563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1249866) q[2];
sx q[2];
rz(-1.5134483) q[2];
sx q[2];
rz(-2.5561257) q[2];
rz(-2.3440907) q[3];
sx q[3];
rz(-1.5944642) q[3];
sx q[3];
rz(-0.40369478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.1095235) q[0];
sx q[0];
rz(-1.7970947) q[0];
sx q[0];
rz(-2.4865785) q[0];
rz(0.0034927448) q[1];
sx q[1];
rz(-0.73492903) q[1];
sx q[1];
rz(3.1266812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92214092) q[0];
sx q[0];
rz(-0.08481124) q[0];
sx q[0];
rz(-0.30973367) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7194233) q[2];
sx q[2];
rz(-0.91277307) q[2];
sx q[2];
rz(-2.7294788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37624371) q[1];
sx q[1];
rz(-2.2606008) q[1];
sx q[1];
rz(2.6756555) q[1];
rz(-pi) q[2];
rz(1.6388325) q[3];
sx q[3];
rz(-2.3792276) q[3];
sx q[3];
rz(-0.12187863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1917176) q[2];
sx q[2];
rz(-1.3363375) q[2];
sx q[2];
rz(0.28676644) q[2];
rz(-0.030335434) q[3];
sx q[3];
rz(-1.560863) q[3];
sx q[3];
rz(-2.0386157) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2716118) q[0];
sx q[0];
rz(-2.5219707) q[0];
sx q[0];
rz(1.9770589) q[0];
rz(-2.8097235) q[1];
sx q[1];
rz(-1.7121366) q[1];
sx q[1];
rz(-1.1480931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6004922) q[0];
sx q[0];
rz(-2.9229188) q[0];
sx q[0];
rz(-1.4906916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0942865) q[2];
sx q[2];
rz(-1.4226573) q[2];
sx q[2];
rz(-2.8830607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0169605) q[1];
sx q[1];
rz(-1.7631044) q[1];
sx q[1];
rz(1.2619727) q[1];
rz(-pi) q[2];
rz(0.049552187) q[3];
sx q[3];
rz(-1.3930071) q[3];
sx q[3];
rz(-2.4018198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0120734) q[2];
sx q[2];
rz(-1.4612863) q[2];
sx q[2];
rz(-2.6626383) q[2];
rz(-2.3524513) q[3];
sx q[3];
rz(-2.456587) q[3];
sx q[3];
rz(-1.3705378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0399465) q[0];
sx q[0];
rz(-2.3532823) q[0];
sx q[0];
rz(2.3729861) q[0];
rz(2.6885314) q[1];
sx q[1];
rz(-2.1078347) q[1];
sx q[1];
rz(2.5624842) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44583089) q[0];
sx q[0];
rz(-0.33798744) q[0];
sx q[0];
rz(0.72622444) q[0];
rz(-pi) q[1];
x q[1];
rz(2.082358) q[2];
sx q[2];
rz(-1.2280012) q[2];
sx q[2];
rz(-2.2276218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.107686) q[1];
sx q[1];
rz(-1.6158069) q[1];
sx q[1];
rz(-2.0652524) q[1];
x q[2];
rz(0.39853951) q[3];
sx q[3];
rz(-1.3966148) q[3];
sx q[3];
rz(-0.70268633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96827489) q[2];
sx q[2];
rz(-1.3710794) q[2];
sx q[2];
rz(-1.0260065) q[2];
rz(-0.46938986) q[3];
sx q[3];
rz(-2.1211076) q[3];
sx q[3];
rz(-0.49218407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7475209) q[0];
sx q[0];
rz(-1.6265656) q[0];
sx q[0];
rz(1.0006022) q[0];
rz(-1.8383149) q[1];
sx q[1];
rz(-2.3293827) q[1];
sx q[1];
rz(-2.2607048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39932809) q[0];
sx q[0];
rz(-2.3652746) q[0];
sx q[0];
rz(2.4180883) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3110169) q[2];
sx q[2];
rz(-2.2143281) q[2];
sx q[2];
rz(2.2011592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92915308) q[1];
sx q[1];
rz(-0.87599659) q[1];
sx q[1];
rz(-1.9966182) q[1];
rz(-0.39805746) q[3];
sx q[3];
rz(-2.0144793) q[3];
sx q[3];
rz(1.8026082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96054968) q[2];
sx q[2];
rz(-2.0281823) q[2];
sx q[2];
rz(-1.95365) q[2];
rz(0.12399593) q[3];
sx q[3];
rz(-1.32722) q[3];
sx q[3];
rz(0.39748642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3666298) q[0];
sx q[0];
rz(-2.9149945) q[0];
sx q[0];
rz(-0.15283787) q[0];
rz(2.9951908) q[1];
sx q[1];
rz(-1.3497137) q[1];
sx q[1];
rz(-2.4375367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.31788) q[0];
sx q[0];
rz(-1.5410265) q[0];
sx q[0];
rz(1.5533621) q[0];
rz(0.91135773) q[2];
sx q[2];
rz(-1.8375085) q[2];
sx q[2];
rz(3.0195723) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2347243) q[1];
sx q[1];
rz(-0.14130302) q[1];
sx q[1];
rz(-1.0060746) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34427963) q[3];
sx q[3];
rz(-2.162291) q[3];
sx q[3];
rz(2.5907093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.081508962) q[2];
sx q[2];
rz(-1.0241822) q[2];
sx q[2];
rz(0.3788968) q[2];
rz(2.1156408) q[3];
sx q[3];
rz(-2.43695) q[3];
sx q[3];
rz(2.0254693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7260471) q[0];
sx q[0];
rz(-1.3274095) q[0];
sx q[0];
rz(2.5782247) q[0];
rz(2.7401961) q[1];
sx q[1];
rz(-2.1386264) q[1];
sx q[1];
rz(-0.67515236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1126661) q[0];
sx q[0];
rz(-0.033551667) q[0];
sx q[0];
rz(1.2262418) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.364351) q[2];
sx q[2];
rz(-2.1690302) q[2];
sx q[2];
rz(-0.75144671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3936972) q[1];
sx q[1];
rz(-1.2505479) q[1];
sx q[1];
rz(-0.25022479) q[1];
x q[2];
rz(-0.37547164) q[3];
sx q[3];
rz(-2.2710147) q[3];
sx q[3];
rz(-0.59798542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99951619) q[2];
sx q[2];
rz(-2.2579305) q[2];
sx q[2];
rz(-0.42259541) q[2];
rz(-2.8790867) q[3];
sx q[3];
rz(-2.0764949) q[3];
sx q[3];
rz(2.0601823) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.066684) q[0];
sx q[0];
rz(-2.1796362) q[0];
sx q[0];
rz(0.99895507) q[0];
rz(-1.2133489) q[1];
sx q[1];
rz(-1.1937001) q[1];
sx q[1];
rz(2.8003069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46743917) q[0];
sx q[0];
rz(-1.7629503) q[0];
sx q[0];
rz(-1.7056607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90088021) q[2];
sx q[2];
rz(-1.7770801) q[2];
sx q[2];
rz(-1.8697949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5707376) q[1];
sx q[1];
rz(-0.26481095) q[1];
sx q[1];
rz(-2.6874506) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.083087) q[3];
sx q[3];
rz(-1.2924177) q[3];
sx q[3];
rz(0.52245058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9294032) q[2];
sx q[2];
rz(-1.5430278) q[2];
sx q[2];
rz(-2.8821778) q[2];
rz(0.2855531) q[3];
sx q[3];
rz(-0.67215896) q[3];
sx q[3];
rz(1.0478421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24006537) q[0];
sx q[0];
rz(-0.49235383) q[0];
sx q[0];
rz(-1.0858076) q[0];
rz(1.0049413) q[1];
sx q[1];
rz(-2.0045547) q[1];
sx q[1];
rz(2.5850632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426612) q[0];
sx q[0];
rz(-0.13971162) q[0];
sx q[0];
rz(1.8675787) q[0];
rz(-2.3031844) q[2];
sx q[2];
rz(-2.072022) q[2];
sx q[2];
rz(-1.7019513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5270015) q[1];
sx q[1];
rz(-1.4785936) q[1];
sx q[1];
rz(0.27669497) q[1];
rz(2.0865165) q[3];
sx q[3];
rz(-2.0175126) q[3];
sx q[3];
rz(1.3019671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13653798) q[2];
sx q[2];
rz(-1.4056861) q[2];
sx q[2];
rz(1.271099) q[2];
rz(0.84730411) q[3];
sx q[3];
rz(-1.3264791) q[3];
sx q[3];
rz(-2.8443851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.51410455) q[0];
sx q[0];
rz(-1.6922373) q[0];
sx q[0];
rz(0.76747146) q[0];
rz(2.1766359) q[1];
sx q[1];
rz(-1.283353) q[1];
sx q[1];
rz(-2.8585785) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3393769) q[0];
sx q[0];
rz(-0.87516038) q[0];
sx q[0];
rz(2.8921769) q[0];
rz(2.709952) q[2];
sx q[2];
rz(-1.5520763) q[2];
sx q[2];
rz(1.4399547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1708888) q[1];
sx q[1];
rz(-1.0701985) q[1];
sx q[1];
rz(-1.1584343) q[1];
x q[2];
rz(-0.55693407) q[3];
sx q[3];
rz(-0.86268988) q[3];
sx q[3];
rz(0.96615893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9877801) q[2];
sx q[2];
rz(-0.38975468) q[2];
sx q[2];
rz(1.7392996) q[2];
rz(0.66579372) q[3];
sx q[3];
rz(-1.8699162) q[3];
sx q[3];
rz(-1.8388892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5182198) q[0];
sx q[0];
rz(-1.8898531) q[0];
sx q[0];
rz(1.3367597) q[0];
rz(1.5917336) q[1];
sx q[1];
rz(-1.5668329) q[1];
sx q[1];
rz(-0.11650539) q[1];
rz(-2.0143853) q[2];
sx q[2];
rz(-2.1090322) q[2];
sx q[2];
rz(2.3354989) q[2];
rz(-2.9339092) q[3];
sx q[3];
rz(-0.80266914) q[3];
sx q[3];
rz(-1.707984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
