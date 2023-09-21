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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5759597) q[0];
sx q[0];
rz(-0.48086777) q[0];
sx q[0];
rz(-0.54310449) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3517411) q[2];
sx q[2];
rz(-2.5918505) q[2];
sx q[2];
rz(-1.654939) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22908902) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(1.2234664) q[1];
x q[2];
rz(0.79521631) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(-2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(2.6057459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8230096) q[0];
sx q[0];
rz(-2.5479925) q[0];
sx q[0];
rz(-1.7943322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8500154) q[2];
sx q[2];
rz(-2.4107286) q[2];
sx q[2];
rz(0.72327327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6260813) q[1];
sx q[1];
rz(-1.1647845) q[1];
sx q[1];
rz(1.8477693) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7137394) q[3];
sx q[3];
rz(-1.5385475) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(0.95834857) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(-0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(-0.025807468) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919825) q[0];
sx q[0];
rz(-1.550807) q[0];
sx q[0];
rz(-1.32094) q[0];
rz(-pi) q[1];
rz(-1.8823207) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(-0.66750079) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4986213) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(0.69795124) q[1];
rz(-0.13173007) q[3];
sx q[3];
rz(-1.9506491) q[3];
sx q[3];
rz(-0.47117885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(0.53612971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49318424) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(2.1302845) q[0];
rz(-pi) q[1];
rz(-2.9085607) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(0.88569966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22524658) q[1];
sx q[1];
rz(-0.9496453) q[1];
sx q[1];
rz(-1.1017373) q[1];
rz(-0.67646497) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(1.3456618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-1.0446576) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(-3.0984745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9263822) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-0.2290639) q[0];
x q[1];
rz(2.8565065) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(1.8269055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.67152126) q[1];
sx q[1];
rz(-2.7579691) q[1];
sx q[1];
rz(0.77312153) q[1];
x q[2];
rz(-0.14857265) q[3];
sx q[3];
rz(-0.79509495) q[3];
sx q[3];
rz(0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(0.56110704) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9583225) q[0];
sx q[0];
rz(-1.1407307) q[0];
sx q[0];
rz(-1.3413315) q[0];
x q[1];
rz(2.6201453) q[2];
sx q[2];
rz(-1.8910742) q[2];
sx q[2];
rz(-2.5495868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42029542) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(-0.23423127) q[1];
x q[2];
rz(0.22535725) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(-2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.51145) q[0];
rz(1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(2.2999433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232988) q[0];
sx q[0];
rz(-2.5677486) q[0];
sx q[0];
rz(0.42563514) q[0];
rz(-pi) q[1];
rz(-0.26288962) q[2];
sx q[2];
rz(-0.6436231) q[2];
sx q[2];
rz(-2.7123244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0824273) q[1];
sx q[1];
rz(-1.6470243) q[1];
sx q[1];
rz(0.04870292) q[1];
x q[2];
rz(0.82960415) q[3];
sx q[3];
rz(-1.7332352) q[3];
sx q[3];
rz(0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324875) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(-1.0743272) q[0];
rz(-pi) q[1];
rz(1.2049963) q[2];
sx q[2];
rz(-1.4636283) q[2];
sx q[2];
rz(-2.0049713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(-2.5820877) q[1];
x q[2];
rz(0.91248625) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(-1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-2.2414482) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8326571) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-0.94690698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98012053) q[0];
sx q[0];
rz(-1.2456018) q[0];
sx q[0];
rz(1.253771) q[0];
rz(-pi) q[1];
rz(2.7526555) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-2.7742085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0526035) q[1];
sx q[1];
rz(-2.7967884) q[1];
sx q[1];
rz(-2.8903972) q[1];
rz(-pi) q[2];
rz(1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(0.63463075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(-2.4865436) q[0];
rz(0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(-2.4972829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5376741) q[0];
sx q[0];
rz(-2.9593421) q[0];
sx q[0];
rz(-0.3084348) q[0];
rz(1.8337433) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(2.3761689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.555968) q[1];
sx q[1];
rz(-1.2871736) q[1];
sx q[1];
rz(-2.7908299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6976922) q[3];
sx q[3];
rz(-1.1392986) q[3];
sx q[3];
rz(2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.775734) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29466378) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-0.12702604) q[2];
sx q[2];
rz(-2.1913678) q[2];
sx q[2];
rz(0.36507228) q[2];
rz(-0.39137822) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
