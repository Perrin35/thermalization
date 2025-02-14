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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(-2.6068249) q[0];
rz(-2.1575902) q[1];
sx q[1];
rz(-0.50552955) q[1];
sx q[1];
rz(-1.2619789) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0875435) q[0];
sx q[0];
rz(-2.4719878) q[0];
sx q[0];
rz(2.7101507) q[0];
rz(-2.1970799) q[2];
sx q[2];
rz(-1.1188095) q[2];
sx q[2];
rz(-2.1747957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3784005) q[1];
sx q[1];
rz(-1.0240558) q[1];
sx q[1];
rz(-1.5101456) q[1];
x q[2];
rz(1.7746968) q[3];
sx q[3];
rz(-1.1887906) q[3];
sx q[3];
rz(2.6873858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24903211) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(1.6585635) q[2];
rz(-0.28462166) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(1.0641789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2738709) q[0];
sx q[0];
rz(-2.4148648) q[0];
sx q[0];
rz(3.100585) q[0];
rz(1.9339336) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(1.4606732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631993) q[0];
sx q[0];
rz(-1.600287) q[0];
sx q[0];
rz(0.14568744) q[0];
x q[1];
rz(2.3744205) q[2];
sx q[2];
rz(-1.1572654) q[2];
sx q[2];
rz(-2.3656379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9750614) q[1];
sx q[1];
rz(-1.6094065) q[1];
sx q[1];
rz(-1.8273942) q[1];
x q[2];
rz(-2.9462141) q[3];
sx q[3];
rz(-2.2535705) q[3];
sx q[3];
rz(-1.186779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4760806) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(-1.7512789) q[2];
rz(3.0604073) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.7056874) q[3];
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
rz(1.3944893) q[0];
sx q[0];
rz(-0.012270027) q[0];
sx q[0];
rz(-2.5315206) q[0];
rz(1.3129781) q[1];
sx q[1];
rz(-0.6956296) q[1];
sx q[1];
rz(2.5909766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3859138) q[0];
sx q[0];
rz(-2.3649594) q[0];
sx q[0];
rz(-2.617111) q[0];
rz(-pi) q[1];
rz(-0.32260334) q[2];
sx q[2];
rz(-1.8292973) q[2];
sx q[2];
rz(2.5532364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50013257) q[1];
sx q[1];
rz(-1.8696897) q[1];
sx q[1];
rz(0.19708339) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9987644) q[3];
sx q[3];
rz(-0.83072829) q[3];
sx q[3];
rz(-0.21634858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47784558) q[2];
sx q[2];
rz(-0.17243324) q[2];
sx q[2];
rz(-2.0752068) q[2];
rz(1.1586698) q[3];
sx q[3];
rz(-1.3830769) q[3];
sx q[3];
rz(0.042044736) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313507) q[0];
sx q[0];
rz(-2.5284335) q[0];
sx q[0];
rz(2.2245275) q[0];
rz(-0.65545583) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(-0.38958946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8389908) q[0];
sx q[0];
rz(-2.9377958) q[0];
sx q[0];
rz(-1.9857282) q[0];
rz(-pi) q[1];
rz(-2.2245313) q[2];
sx q[2];
rz(-2.3417763) q[2];
sx q[2];
rz(0.66674495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45976156) q[1];
sx q[1];
rz(-1.1696522) q[1];
sx q[1];
rz(-0.031706867) q[1];
rz(-2.5960931) q[3];
sx q[3];
rz(-2.0974468) q[3];
sx q[3];
rz(2.5712476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3098844) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(1.3523098) q[2];
rz(2.3934707) q[3];
sx q[3];
rz(-1.8375405) q[3];
sx q[3];
rz(1.6149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-0.56981531) q[0];
sx q[0];
rz(-1.3414398) q[0];
rz(-2.1603284) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(-0.70995465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34559743) q[0];
sx q[0];
rz(-1.5516722) q[0];
sx q[0];
rz(-0.32581331) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1469028) q[2];
sx q[2];
rz(-1.369595) q[2];
sx q[2];
rz(-0.76801571) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4932351) q[1];
sx q[1];
rz(-1.7410189) q[1];
sx q[1];
rz(1.3551718) q[1];
x q[2];
rz(2.5437935) q[3];
sx q[3];
rz(-1.6493268) q[3];
sx q[3];
rz(-0.30818916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-2.3812582) q[2];
sx q[2];
rz(2.23488) q[2];
rz(-0.18236154) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(-2.2831634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3458503) q[0];
sx q[0];
rz(-0.95473552) q[0];
sx q[0];
rz(2.9398651) q[0];
rz(-1.7851104) q[1];
sx q[1];
rz(-0.67142612) q[1];
sx q[1];
rz(1.9904402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042554458) q[0];
sx q[0];
rz(-1.3622314) q[0];
sx q[0];
rz(2.0046356) q[0];
rz(-pi) q[1];
rz(0.62080748) q[2];
sx q[2];
rz(-0.70638958) q[2];
sx q[2];
rz(2.8548129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5804739) q[1];
sx q[1];
rz(-2.4137133) q[1];
sx q[1];
rz(2.5125458) q[1];
rz(-pi) q[2];
rz(-2.7490691) q[3];
sx q[3];
rz(-0.54504921) q[3];
sx q[3];
rz(-0.98859331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1543033) q[2];
sx q[2];
rz(-2.5550948) q[2];
sx q[2];
rz(0.34745535) q[2];
rz(2.3197428) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(-1.6177026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007941) q[0];
sx q[0];
rz(-2.4947385) q[0];
sx q[0];
rz(-1.1232173) q[0];
rz(2.7128291) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(2.9926328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4377047) q[0];
sx q[0];
rz(-2.0876679) q[0];
sx q[0];
rz(0.025512841) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8455714) q[2];
sx q[2];
rz(-1.3775577) q[2];
sx q[2];
rz(-3.0770965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91576112) q[1];
sx q[1];
rz(-1.718545) q[1];
sx q[1];
rz(-0.76051449) q[1];
rz(1.1650208) q[3];
sx q[3];
rz(-1.6029442) q[3];
sx q[3];
rz(-0.11628499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.011977) q[2];
sx q[2];
rz(-2.7836383) q[2];
sx q[2];
rz(0.32148662) q[2];
rz(2.0251515) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(-1.4139676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1534934) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(-0.022911428) q[0];
rz(-1.5921536) q[1];
sx q[1];
rz(-2.0982845) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003636) q[0];
sx q[0];
rz(-2.1080906) q[0];
sx q[0];
rz(-2.8365718) q[0];
rz(1.171807) q[2];
sx q[2];
rz(-2.3697457) q[2];
sx q[2];
rz(-3.0400624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0523394) q[1];
sx q[1];
rz(-1.2683378) q[1];
sx q[1];
rz(-0.41771981) q[1];
rz(-0.096209196) q[3];
sx q[3];
rz(-2.2819977) q[3];
sx q[3];
rz(-0.32207738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93078485) q[2];
sx q[2];
rz(-2.942694) q[2];
sx q[2];
rz(0.99736324) q[2];
rz(1.2584244) q[3];
sx q[3];
rz(-1.9505898) q[3];
sx q[3];
rz(1.9112126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9835984) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(-0.47437814) q[0];
rz(0.55167088) q[1];
sx q[1];
rz(-0.76849476) q[1];
sx q[1];
rz(-0.65753585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6303355) q[0];
sx q[0];
rz(-1.5310107) q[0];
sx q[0];
rz(-1.5471094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9766519) q[2];
sx q[2];
rz(-1.4144344) q[2];
sx q[2];
rz(2.641822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8213013) q[1];
sx q[1];
rz(-1.5870915) q[1];
sx q[1];
rz(1.2024131) q[1];
rz(3.0364584) q[3];
sx q[3];
rz(-1.7145559) q[3];
sx q[3];
rz(-1.2105699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9115596) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(1.3113021) q[2];
rz(0.71410549) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(-0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(-2.5332992) q[0];
rz(-0.81781203) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(0.83293319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8084106) q[0];
sx q[0];
rz(-1.6224553) q[0];
sx q[0];
rz(2.7201061) q[0];
rz(2.8274546) q[2];
sx q[2];
rz(-1.8153662) q[2];
sx q[2];
rz(0.42942522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45730879) q[1];
sx q[1];
rz(-0.59694081) q[1];
sx q[1];
rz(0.2376016) q[1];
rz(-pi) q[2];
rz(1.5394443) q[3];
sx q[3];
rz(-0.97797457) q[3];
sx q[3];
rz(-1.5167936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9587162) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(0.58615169) q[2];
rz(0.90306774) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(0.085845145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387954) q[0];
sx q[0];
rz(-2.0113404) q[0];
sx q[0];
rz(2.197862) q[0];
rz(2.0441652) q[1];
sx q[1];
rz(-2.8891017) q[1];
sx q[1];
rz(-0.15695922) q[1];
rz(-1.912276) q[2];
sx q[2];
rz(-1.7021644) q[2];
sx q[2];
rz(1.1973039) q[2];
rz(1.0020574) q[3];
sx q[3];
rz(-2.3593223) q[3];
sx q[3];
rz(0.13025688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
