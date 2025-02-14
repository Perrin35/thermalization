OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.77312624) q[0];
sx q[0];
rz(2.3946895) q[0];
sx q[0];
rz(11.725732) q[0];
rz(3.2631915) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(12.267332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6071607) q[0];
sx q[0];
rz(-2.2272155) q[0];
sx q[0];
rz(0.54404152) q[0];
rz(-pi) q[1];
rz(-1.0844803) q[2];
sx q[2];
rz(-1.5019226) q[2];
sx q[2];
rz(1.8191847) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2192397) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(2.8994865) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2848008) q[3];
sx q[3];
rz(-2.6702849) q[3];
sx q[3];
rz(2.907674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3806939) q[2];
sx q[2];
rz(-1.0502522) q[2];
sx q[2];
rz(0.32763457) q[2];
rz(-1.3753752) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.37680092) q[0];
sx q[0];
rz(-1.5400274) q[0];
sx q[0];
rz(-0.85533992) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(-0.45113742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6159003) q[0];
sx q[0];
rz(-0.85291686) q[0];
sx q[0];
rz(2.9327716) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7629659) q[2];
sx q[2];
rz(-0.8079257) q[2];
sx q[2];
rz(2.5469123) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6169713) q[1];
sx q[1];
rz(-2.0071623) q[1];
sx q[1];
rz(-1.7625767) q[1];
x q[2];
rz(-1.9782009) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(0.16678424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96192876) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(-0.12953225) q[2];
rz(-2.9642963) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(-3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(1.7496846) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-0.55927292) q[0];
rz(3.0468805) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(2.6643378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9608979) q[0];
sx q[0];
rz(-1.8261693) q[0];
sx q[0];
rz(-1.9810505) q[0];
x q[1];
rz(-2.7626286) q[2];
sx q[2];
rz(-2.7550089) q[2];
sx q[2];
rz(1.6340337) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92465688) q[1];
sx q[1];
rz(-2.3291596) q[1];
sx q[1];
rz(-1.8820018) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84945143) q[3];
sx q[3];
rz(-1.9799383) q[3];
sx q[3];
rz(-2.9144998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7462848) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(1.2578472) q[3];
sx q[3];
rz(-1.427429) q[3];
sx q[3];
rz(-3.0959082) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78087085) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(-1.038653) q[0];
rz(-2.5301798) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(2.2183653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99902376) q[0];
sx q[0];
rz(-1.2272738) q[0];
sx q[0];
rz(-0.96238636) q[0];
rz(-pi) q[1];
rz(1.5973041) q[2];
sx q[2];
rz(-2.7237281) q[2];
sx q[2];
rz(-0.94965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9070396) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(-1.4437463) q[1];
rz(-pi) q[2];
rz(-1.3096894) q[3];
sx q[3];
rz(-0.93702836) q[3];
sx q[3];
rz(1.5719617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63127798) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-2.5941217) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(-2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15544686) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(0.11988457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6836707) q[0];
sx q[0];
rz(-0.73212762) q[0];
sx q[0];
rz(-1.2326272) q[0];
rz(-2.2861459) q[2];
sx q[2];
rz(-2.3513633) q[2];
sx q[2];
rz(1.9270735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6189055) q[1];
sx q[1];
rz(-0.22929103) q[1];
sx q[1];
rz(0.35405901) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0644887) q[3];
sx q[3];
rz(-2.1112006) q[3];
sx q[3];
rz(-0.04549724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1025866) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(-2.8809123) q[2];
rz(-0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(-2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6108625) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(2.196905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8588036) q[0];
sx q[0];
rz(-1.7886046) q[0];
sx q[0];
rz(-2.5522638) q[0];
rz(-pi) q[1];
rz(-0.86017365) q[2];
sx q[2];
rz(-0.80931907) q[2];
sx q[2];
rz(0.10242538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78996736) q[1];
sx q[1];
rz(-1.891544) q[1];
sx q[1];
rz(1.411639) q[1];
rz(-pi) q[2];
rz(1.764749) q[3];
sx q[3];
rz(-2.5911281) q[3];
sx q[3];
rz(0.46904072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(0.043924335) q[2];
rz(1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.6759492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2783022) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(2.6712766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4620275) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(-0.26058773) q[0];
x q[1];
rz(1.1675179) q[2];
sx q[2];
rz(-1.2263311) q[2];
sx q[2];
rz(-3.0967876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1086169) q[1];
sx q[1];
rz(-1.2101296) q[1];
sx q[1];
rz(-0.74957871) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10736088) q[3];
sx q[3];
rz(-0.11365644) q[3];
sx q[3];
rz(-1.3945127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(0.31141591) q[2];
rz(-0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(-0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-2.2054963) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5514126) q[0];
sx q[0];
rz(-1.0820933) q[0];
sx q[0];
rz(1.9103229) q[0];
rz(-pi) q[1];
rz(-3.1390433) q[2];
sx q[2];
rz(-1.7685206) q[2];
sx q[2];
rz(-1.6127123) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0691904) q[1];
sx q[1];
rz(-1.6399929) q[1];
sx q[1];
rz(-1.9131768) q[1];
rz(-pi) q[2];
rz(-1.3948729) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(-1.3678838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(1.4754971) q[2];
rz(2.3462319) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(-3.0812145) q[0];
rz(-2.9810442) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(0.9789595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5125324) q[0];
sx q[0];
rz(-0.94659014) q[0];
sx q[0];
rz(1.4275622) q[0];
rz(-pi) q[1];
rz(-1.2376182) q[2];
sx q[2];
rz(-2.4314483) q[2];
sx q[2];
rz(1.4613446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97040788) q[1];
sx q[1];
rz(-1.1029585) q[1];
sx q[1];
rz(2.4572608) q[1];
rz(-pi) q[2];
x q[2];
rz(0.03421182) q[3];
sx q[3];
rz(-2.1497823) q[3];
sx q[3];
rz(-1.2466696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0673361) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(2.5439751) q[3];
sx q[3];
rz(-1.7643192) q[3];
sx q[3];
rz(-0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.446796) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(0.52892518) q[0];
rz(2.9534598) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(0.58427748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54160252) q[0];
sx q[0];
rz(-1.3253821) q[0];
sx q[0];
rz(2.8757069) q[0];
rz(-pi) q[1];
rz(-1.9599914) q[2];
sx q[2];
rz(-1.1616594) q[2];
sx q[2];
rz(-1.429806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53946282) q[1];
sx q[1];
rz(-1.6136843) q[1];
sx q[1];
rz(2.7320288) q[1];
rz(-0.27354555) q[3];
sx q[3];
rz(-2.1714032) q[3];
sx q[3];
rz(2.724444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75954413) q[2];
sx q[2];
rz(-1.0039165) q[2];
sx q[2];
rz(3.0697401) q[2];
rz(2.1182649) q[3];
sx q[3];
rz(-1.5724678) q[3];
sx q[3];
rz(-0.48880997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857984) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(-0.51495348) q[2];
sx q[2];
rz(-1.6390159) q[2];
sx q[2];
rz(-1.9255571) q[2];
rz(0.16719462) q[3];
sx q[3];
rz(-1.7291369) q[3];
sx q[3];
rz(-0.72469934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
