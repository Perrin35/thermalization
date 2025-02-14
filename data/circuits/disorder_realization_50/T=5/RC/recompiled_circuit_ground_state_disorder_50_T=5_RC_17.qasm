OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9775554) q[0];
sx q[0];
rz(-1.2241751) q[0];
sx q[0];
rz(1.2807711) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(-2.8746936) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3506265) q[0];
sx q[0];
rz(-1.4746951) q[0];
sx q[0];
rz(-0.093856974) q[0];
rz(-pi) q[1];
rz(1.880015) q[2];
sx q[2];
rz(-0.71239939) q[2];
sx q[2];
rz(1.0605896) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34632979) q[1];
sx q[1];
rz(-0.85546934) q[1];
sx q[1];
rz(-0.99154179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6736753) q[3];
sx q[3];
rz(-1.130338) q[3];
sx q[3];
rz(-2.6552424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12307564) q[2];
sx q[2];
rz(-0.96394959) q[2];
sx q[2];
rz(0.71259552) q[2];
rz(0.16768843) q[3];
sx q[3];
rz(-2.2919787) q[3];
sx q[3];
rz(0.4445506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0406822) q[0];
sx q[0];
rz(-1.7829144) q[0];
sx q[0];
rz(-1.4605301) q[0];
rz(-1.679861) q[1];
sx q[1];
rz(-2.0748383) q[1];
sx q[1];
rz(1.1163968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8212962) q[0];
sx q[0];
rz(-1.5954994) q[0];
sx q[0];
rz(-1.5276093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15123151) q[2];
sx q[2];
rz(-1.0200715) q[2];
sx q[2];
rz(-0.085892396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18552407) q[1];
sx q[1];
rz(-1.1949693) q[1];
sx q[1];
rz(1.2898499) q[1];
rz(1.7061911) q[3];
sx q[3];
rz(-1.1639126) q[3];
sx q[3];
rz(-0.02449206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1560893) q[2];
sx q[2];
rz(-0.49571338) q[2];
sx q[2];
rz(2.8893341) q[2];
rz(-1.0559399) q[3];
sx q[3];
rz(-2.2838433) q[3];
sx q[3];
rz(1.0004388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58241874) q[0];
sx q[0];
rz(-2.2370339) q[0];
sx q[0];
rz(0.82758033) q[0];
rz(-0.52455348) q[1];
sx q[1];
rz(-1.2453715) q[1];
sx q[1];
rz(-2.1771199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37999718) q[0];
sx q[0];
rz(-2.5598713) q[0];
sx q[0];
rz(1.6935478) q[0];
rz(-1.9461412) q[2];
sx q[2];
rz(-1.4895194) q[2];
sx q[2];
rz(0.92513212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4276994) q[1];
sx q[1];
rz(-1.6203708) q[1];
sx q[1];
rz(-0.63062856) q[1];
rz(-pi) q[2];
rz(3.019624) q[3];
sx q[3];
rz(-0.44443529) q[3];
sx q[3];
rz(-1.1201064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39782897) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(1.1986097) q[2];
rz(1.4926636) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(-2.5627513) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156133) q[0];
sx q[0];
rz(-0.15376832) q[0];
sx q[0];
rz(-1.2157259) q[0];
rz(0.99023306) q[1];
sx q[1];
rz(-0.74140397) q[1];
sx q[1];
rz(-0.56891099) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.194806) q[0];
sx q[0];
rz(-1.0005568) q[0];
sx q[0];
rz(-2.8034745) q[0];
rz(-pi) q[1];
rz(-2.8773035) q[2];
sx q[2];
rz(-1.2131872) q[2];
sx q[2];
rz(-0.25928318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7987631) q[1];
sx q[1];
rz(-1.1278648) q[1];
sx q[1];
rz(2.8196067) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91671555) q[3];
sx q[3];
rz(-0.93610969) q[3];
sx q[3];
rz(-3.0460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8625921) q[2];
sx q[2];
rz(-0.96070868) q[2];
sx q[2];
rz(-0.29435364) q[2];
rz(-2.477008) q[3];
sx q[3];
rz(-0.75988257) q[3];
sx q[3];
rz(-1.6871281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45254529) q[0];
sx q[0];
rz(-1.713546) q[0];
sx q[0];
rz(2.8277165) q[0];
rz(0.9017871) q[1];
sx q[1];
rz(-1.3548464) q[1];
sx q[1];
rz(-1.4656167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70756992) q[0];
sx q[0];
rz(-1.5450053) q[0];
sx q[0];
rz(-1.6898193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9889262) q[2];
sx q[2];
rz(-0.2731495) q[2];
sx q[2];
rz(1.404431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5295388) q[1];
sx q[1];
rz(-1.1837237) q[1];
sx q[1];
rz(-1.6516634) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63057804) q[3];
sx q[3];
rz(-1.6007716) q[3];
sx q[3];
rz(1.1512299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9471211) q[2];
sx q[2];
rz(-0.7496382) q[2];
sx q[2];
rz(2.0571902) q[2];
rz(0.32367745) q[3];
sx q[3];
rz(-0.71666986) q[3];
sx q[3];
rz(-2.9173541) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502458) q[0];
sx q[0];
rz(-2.8960189) q[0];
sx q[0];
rz(0.41159758) q[0];
rz(0.72548524) q[1];
sx q[1];
rz(-1.2440224) q[1];
sx q[1];
rz(1.4010319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66564225) q[0];
sx q[0];
rz(-0.67806292) q[0];
sx q[0];
rz(-0.31291385) q[0];
rz(-pi) q[1];
rz(-2.4011924) q[2];
sx q[2];
rz(-1.6639189) q[2];
sx q[2];
rz(1.4880866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2620777) q[1];
sx q[1];
rz(-1.8873155) q[1];
sx q[1];
rz(0.36088093) q[1];
rz(-pi) q[2];
rz(0.45127604) q[3];
sx q[3];
rz(-1.5022931) q[3];
sx q[3];
rz(0.20718224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97658849) q[2];
sx q[2];
rz(-2.5133331) q[2];
sx q[2];
rz(-1.4264301) q[2];
rz(0.12380883) q[3];
sx q[3];
rz(-2.0721469) q[3];
sx q[3];
rz(0.4933221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2774778) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(-1.3364828) q[0];
rz(-2.2026964) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(0.28405651) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2920609) q[0];
sx q[0];
rz(-0.96905989) q[0];
sx q[0];
rz(1.1547778) q[0];
x q[1];
rz(2.9781173) q[2];
sx q[2];
rz(-1.3473841) q[2];
sx q[2];
rz(-0.3439518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5732481) q[1];
sx q[1];
rz(-1.580984) q[1];
sx q[1];
rz(1.5393901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18064336) q[3];
sx q[3];
rz(-1.4339281) q[3];
sx q[3];
rz(-1.6302072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17860086) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(2.2006688) q[2];
rz(0.00087794463) q[3];
sx q[3];
rz(-1.9160756) q[3];
sx q[3];
rz(-3.0555449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451999) q[0];
sx q[0];
rz(-0.89458507) q[0];
sx q[0];
rz(0.17564242) q[0];
rz(1.2200217) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(2.2697935) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0891465) q[0];
sx q[0];
rz(-1.8348357) q[0];
sx q[0];
rz(-2.9857078) q[0];
x q[1];
rz(-0.37105889) q[2];
sx q[2];
rz(-2.3239229) q[2];
sx q[2];
rz(2.121986) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6260274) q[1];
sx q[1];
rz(-0.62735081) q[1];
sx q[1];
rz(-2.8684142) q[1];
rz(-pi) q[2];
rz(1.4967887) q[3];
sx q[3];
rz(-1.0153061) q[3];
sx q[3];
rz(0.018270103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4153727) q[2];
sx q[2];
rz(-1.2803593) q[2];
sx q[2];
rz(1.0225164) q[2];
rz(0.38698777) q[3];
sx q[3];
rz(-1.3849247) q[3];
sx q[3];
rz(-2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047394) q[0];
sx q[0];
rz(-1.1664718) q[0];
sx q[0];
rz(2.8785896) q[0];
rz(-2.8291342) q[1];
sx q[1];
rz(-1.0954906) q[1];
sx q[1];
rz(1.9591263) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179494) q[0];
sx q[0];
rz(-1.765476) q[0];
sx q[0];
rz(-1.0561835) q[0];
rz(-2.1698007) q[2];
sx q[2];
rz(-2.3070719) q[2];
sx q[2];
rz(2.4311299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14543372) q[1];
sx q[1];
rz(-2.7163032) q[1];
sx q[1];
rz(-0.49844679) q[1];
rz(-pi) q[2];
rz(-0.5236756) q[3];
sx q[3];
rz(-2.0770235) q[3];
sx q[3];
rz(2.5158213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.075228127) q[2];
sx q[2];
rz(-1.3508537) q[2];
sx q[2];
rz(2.8933375) q[2];
rz(-2.9634641) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(-2.3586912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42661509) q[0];
sx q[0];
rz(-2.7474032) q[0];
sx q[0];
rz(-2.07975) q[0];
rz(-1.8408403) q[1];
sx q[1];
rz(-1.6164833) q[1];
sx q[1];
rz(-1.1311857) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0164364) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(2.6317276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50297351) q[2];
sx q[2];
rz(-1.1736693) q[2];
sx q[2];
rz(-3.0850038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6415875) q[1];
sx q[1];
rz(-1.4785071) q[1];
sx q[1];
rz(-2.7422043) q[1];
rz(-2.5401073) q[3];
sx q[3];
rz(-1.7643398) q[3];
sx q[3];
rz(-2.7880965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0807557) q[2];
sx q[2];
rz(-2.8751825) q[2];
sx q[2];
rz(2.5401435) q[2];
rz(-2.1001749) q[3];
sx q[3];
rz(-1.4789378) q[3];
sx q[3];
rz(1.3781594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9926485) q[0];
sx q[0];
rz(-2.5288378) q[0];
sx q[0];
rz(0.80768325) q[0];
rz(0.07769892) q[1];
sx q[1];
rz(-0.54010375) q[1];
sx q[1];
rz(1.7355951) q[1];
rz(2.7355657) q[2];
sx q[2];
rz(-2.1320504) q[2];
sx q[2];
rz(-0.8036094) q[2];
rz(-2.349825) q[3];
sx q[3];
rz(-2.0590012) q[3];
sx q[3];
rz(1.7760359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
