OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.565633) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(0.54310449) q[0];
rz(-pi) q[1];
rz(1.3517411) q[2];
sx q[2];
rz(-2.5918505) q[2];
sx q[2];
rz(-1.4866536) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5696722) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(0.7260679) q[1];
rz(-pi) q[2];
rz(-0.79521631) q[3];
sx q[3];
rz(-1.8024369) q[3];
sx q[3];
rz(0.29602805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(1.583064) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-2.6057459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5553404) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(2.9931086) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8230121) q[2];
sx q[2];
rz(-0.87715845) q[2];
sx q[2];
rz(2.0351621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51551137) q[1];
sx q[1];
rz(-1.1647845) q[1];
sx q[1];
rz(1.8477693) q[1];
x q[2];
rz(1.7934947) q[3];
sx q[3];
rz(-0.14651146) q[3];
sx q[3];
rz(3.1080064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(2.9911175) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31608554) q[0];
sx q[0];
rz(-1.320991) q[0];
sx q[0];
rz(-0.020629701) q[0];
rz(-1.843156) q[2];
sx q[2];
rz(-2.8189427) q[2];
sx q[2];
rz(-1.9793561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7681231) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(2.5440689) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8886391) q[3];
sx q[3];
rz(-0.40099537) q[3];
sx q[3];
rz(-2.3272115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-0.53612971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25585184) q[0];
sx q[0];
rz(-2.2478074) q[0];
sx q[0];
rz(0.52744249) q[0];
rz(-0.23303194) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(2.255893) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.199898) q[1];
sx q[1];
rz(-2.3823793) q[1];
sx q[1];
rz(-2.5782176) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36066182) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(0.043118127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9163141) q[0];
sx q[0];
rz(-1.3818704) q[0];
sx q[0];
rz(0.96154763) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28508614) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(1.8269055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16380285) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(0.28122854) q[1];
x q[2];
rz(2.99302) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9583225) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(-1.8002611) q[0];
rz(-pi) q[1];
rz(-2.6201453) q[2];
sx q[2];
rz(-1.8910742) q[2];
sx q[2];
rz(2.5495868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7212972) q[1];
sx q[1];
rz(-1.2696206) q[1];
sx q[1];
rz(-2.9073614) q[1];
rz(-1.7644464) q[3];
sx q[3];
rz(-0.87246694) q[3];
sx q[3];
rz(0.26102548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(-0.84164936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8182939) q[0];
sx q[0];
rz(-2.5677486) q[0];
sx q[0];
rz(-0.42563514) q[0];
rz(-pi) q[1];
rz(-1.3782578) q[2];
sx q[2];
rz(-0.95270573) q[2];
sx q[2];
rz(-2.3877909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51534286) q[1];
sx q[1];
rz(-1.5222349) q[1];
sx q[1];
rz(1.6471144) q[1];
rz(-0.21861403) q[3];
sx q[3];
rz(-0.84158763) q[3];
sx q[3];
rz(1.5515755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(2.7239674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091052) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(2.0672654) q[0];
rz(-pi) q[1];
rz(-0.11469658) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(-2.7483658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(-2.5820877) q[1];
rz(-0.29406677) q[3];
sx q[3];
rz(-1.2122279) q[3];
sx q[3];
rz(-2.8420574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(-0.67388326) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(0.62966627) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(2.1946857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6553584) q[0];
sx q[0];
rz(-1.2709193) q[0];
sx q[0];
rz(0.34098682) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7526555) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-2.7742085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0526035) q[1];
sx q[1];
rz(-2.7967884) q[1];
sx q[1];
rz(-0.25119541) q[1];
rz(-pi) q[2];
rz(-1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(-0.63463075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(-2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-0.64430976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6039186) q[0];
sx q[0];
rz(-2.9593421) q[0];
sx q[0];
rz(-0.3084348) q[0];
rz(-pi) q[1];
rz(-0.23004736) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(-2.7237797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2584553) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(1.2698445) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82018567) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(-1.5370777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(2.541686) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(1.3953801) q[2];
sx q[2];
rz(-0.63175628) q[2];
sx q[2];
rz(0.14887688) q[2];
rz(-2.0471845) q[3];
sx q[3];
rz(-0.73275685) q[3];
sx q[3];
rz(-2.3263596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];