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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(0.59803522) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(2.476517) q[1];
sx q[1];
rz(9.2718931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0273697) q[0];
sx q[0];
rz(-1.4589757) q[0];
sx q[0];
rz(0.080145135) q[0];
rz(-pi) q[1];
x q[1];
rz(2.653605) q[2];
sx q[2];
rz(-1.5190912) q[2];
sx q[2];
rz(-2.0590797) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8522332) q[1];
sx q[1];
rz(-0.24670641) q[1];
sx q[1];
rz(1.0323332) q[1];
rz(-pi) q[2];
rz(1.5737719) q[3];
sx q[3];
rz(-1.5248393) q[3];
sx q[3];
rz(2.9480235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87024706) q[2];
sx q[2];
rz(-0.18834867) q[2];
sx q[2];
rz(1.9625473) q[2];
rz(-2.7315268) q[3];
sx q[3];
rz(-2.0867917) q[3];
sx q[3];
rz(2.2470391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9533185) q[0];
sx q[0];
rz(-0.66903791) q[0];
sx q[0];
rz(2.8642995) q[0];
rz(-2.9412728) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(-3.0576113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7070933) q[0];
sx q[0];
rz(-1.3719808) q[0];
sx q[0];
rz(-1.0931703) q[0];
x q[1];
rz(1.7189156) q[2];
sx q[2];
rz(-2.2930055) q[2];
sx q[2];
rz(-0.87004694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9878432) q[1];
sx q[1];
rz(-1.1293518) q[1];
sx q[1];
rz(0.50905871) q[1];
rz(-1.8324319) q[3];
sx q[3];
rz(-1.0519093) q[3];
sx q[3];
rz(2.6760657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1818485) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(0.7461156) q[2];
rz(-0.6212081) q[3];
sx q[3];
rz(-2.4977903) q[3];
sx q[3];
rz(1.3889036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32034945) q[0];
sx q[0];
rz(-2.4140883) q[0];
sx q[0];
rz(1.9387091) q[0];
rz(-0.41162833) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(-1.113755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6516371) q[0];
sx q[0];
rz(-2.2764766) q[0];
sx q[0];
rz(-0.56244992) q[0];
rz(-0.97881563) q[2];
sx q[2];
rz(-2.4360058) q[2];
sx q[2];
rz(-1.6340812) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74677709) q[1];
sx q[1];
rz(-1.703843) q[1];
sx q[1];
rz(2.2486671) q[1];
x q[2];
rz(-0.57856929) q[3];
sx q[3];
rz(-0.95953836) q[3];
sx q[3];
rz(1.3276154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0456298) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(0.80471936) q[3];
sx q[3];
rz(-1.8685124) q[3];
sx q[3];
rz(-1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74748892) q[0];
sx q[0];
rz(-1.8590314) q[0];
sx q[0];
rz(1.6266741) q[0];
rz(-2.6331242) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(1.633684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8430458) q[0];
sx q[0];
rz(-2.9905479) q[0];
sx q[0];
rz(-1.5848978) q[0];
x q[1];
rz(3.0005387) q[2];
sx q[2];
rz(-1.3562437) q[2];
sx q[2];
rz(0.44520865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8123618) q[1];
sx q[1];
rz(-0.32608247) q[1];
sx q[1];
rz(1.2124168) q[1];
rz(-0.42801492) q[3];
sx q[3];
rz(-0.95967442) q[3];
sx q[3];
rz(-2.1111272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7278829) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(-0.74753648) q[2];
rz(1.9109292) q[3];
sx q[3];
rz(-2.3321798) q[3];
sx q[3];
rz(-0.49720732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43664765) q[0];
sx q[0];
rz(-1.1282938) q[0];
sx q[0];
rz(1.4592015) q[0];
rz(2.7875426) q[1];
sx q[1];
rz(-0.67146462) q[1];
sx q[1];
rz(-0.57463247) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1099691) q[0];
sx q[0];
rz(-1.4935483) q[0];
sx q[0];
rz(-1.5573274) q[0];
x q[1];
rz(-2.3726497) q[2];
sx q[2];
rz(-1.7684325) q[2];
sx q[2];
rz(2.1773424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71855356) q[1];
sx q[1];
rz(-1.3438112) q[1];
sx q[1];
rz(1.6237698) q[1];
x q[2];
rz(1.7754885) q[3];
sx q[3];
rz(-1.6833971) q[3];
sx q[3];
rz(-1.7623869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5057126) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(-0.20225784) q[2];
rz(-2.108719) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(-3.0752944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.031484) q[0];
sx q[0];
rz(-2.7543572) q[0];
sx q[0];
rz(-2.7701344) q[0];
rz(0.29608852) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(2.9782226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0609157) q[0];
sx q[0];
rz(-1.0759584) q[0];
sx q[0];
rz(-1.3617152) q[0];
rz(2.4642015) q[2];
sx q[2];
rz(-2.4755726) q[2];
sx q[2];
rz(0.99413727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5529873) q[1];
sx q[1];
rz(-1.3901911) q[1];
sx q[1];
rz(-0.44579472) q[1];
x q[2];
rz(0.87518163) q[3];
sx q[3];
rz(-1.6871678) q[3];
sx q[3];
rz(2.5494408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0614193) q[2];
sx q[2];
rz(-0.27013186) q[2];
sx q[2];
rz(-2.488625) q[2];
rz(2.669892) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(-2.4145943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93040526) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(-2.1042673) q[0];
rz(0.51517454) q[1];
sx q[1];
rz(-1.8637135) q[1];
sx q[1];
rz(-1.3264664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4194834) q[0];
sx q[0];
rz(-2.4631755) q[0];
sx q[0];
rz(2.2161311) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7896129) q[2];
sx q[2];
rz(-1.6087785) q[2];
sx q[2];
rz(-0.69881907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3390914) q[1];
sx q[1];
rz(-2.3330124) q[1];
sx q[1];
rz(-1.436961) q[1];
x q[2];
rz(0.063422261) q[3];
sx q[3];
rz(-1.2655971) q[3];
sx q[3];
rz(2.2535107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6933763) q[2];
sx q[2];
rz(-0.46458149) q[2];
sx q[2];
rz(-0.32619897) q[2];
rz(2.163573) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(-0.090156468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4433032) q[0];
sx q[0];
rz(-2.8841618) q[0];
sx q[0];
rz(0.28939104) q[0];
rz(3.1293213) q[1];
sx q[1];
rz(-0.27996501) q[1];
sx q[1];
rz(1.6562921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78083437) q[0];
sx q[0];
rz(-0.88874431) q[0];
sx q[0];
rz(-2.2498796) q[0];
x q[1];
rz(-0.22356914) q[2];
sx q[2];
rz(-1.5180598) q[2];
sx q[2];
rz(1.2745672) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43013182) q[1];
sx q[1];
rz(-1.7727114) q[1];
sx q[1];
rz(1.422775) q[1];
rz(-0.9910219) q[3];
sx q[3];
rz(-0.54383792) q[3];
sx q[3];
rz(0.71292927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54899186) q[2];
sx q[2];
rz(-1.5697378) q[2];
sx q[2];
rz(0.17295095) q[2];
rz(-1.3744099) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(1.4779199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7462815) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(0.81004274) q[0];
rz(-2.7250302) q[1];
sx q[1];
rz(-2.3655393) q[1];
sx q[1];
rz(2.7966444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0121794) q[0];
sx q[0];
rz(-2.4218771) q[0];
sx q[0];
rz(0.1869633) q[0];
x q[1];
rz(1.0428997) q[2];
sx q[2];
rz(-0.32020346) q[2];
sx q[2];
rz(-2.7803583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7564063) q[1];
sx q[1];
rz(-1.1118206) q[1];
sx q[1];
rz(2.7800757) q[1];
x q[2];
rz(-1.6900605) q[3];
sx q[3];
rz(-1.4550617) q[3];
sx q[3];
rz(-1.3758462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42744669) q[2];
sx q[2];
rz(-0.88627187) q[2];
sx q[2];
rz(0.084058849) q[2];
rz(0.90703026) q[3];
sx q[3];
rz(-1.4182988) q[3];
sx q[3];
rz(-2.1840054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3323988) q[0];
sx q[0];
rz(-3.0539303) q[0];
sx q[0];
rz(2.8293389) q[0];
rz(0.23278438) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(1.2871453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96583073) q[0];
sx q[0];
rz(-1.2487808) q[0];
sx q[0];
rz(-2.470592) q[0];
x q[1];
rz(-3.078545) q[2];
sx q[2];
rz(-2.5015066) q[2];
sx q[2];
rz(0.62713059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47984581) q[1];
sx q[1];
rz(-1.8150618) q[1];
sx q[1];
rz(2.0354009) q[1];
x q[2];
rz(-0.62508836) q[3];
sx q[3];
rz(-1.2660789) q[3];
sx q[3];
rz(2.8947322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1708258) q[2];
sx q[2];
rz(-1.2831186) q[2];
sx q[2];
rz(-1.4975366) q[2];
rz(2.0434642) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(-0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.004414) q[0];
sx q[0];
rz(-1.7902086) q[0];
sx q[0];
rz(-0.93850346) q[0];
rz(-1.3068403) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(-2.2597792) q[2];
sx q[2];
rz(-1.6826993) q[2];
sx q[2];
rz(-0.8347389) q[2];
rz(-2.2372237) q[3];
sx q[3];
rz(-1.4587797) q[3];
sx q[3];
rz(0.070894632) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
