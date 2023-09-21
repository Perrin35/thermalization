OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(4.1783279) q[1];
sx q[1];
rz(8.7611603) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27204313) q[0];
sx q[0];
rz(-0.46177319) q[0];
sx q[0];
rz(-0.053134993) q[0];
rz(-pi) q[1];
rz(0.3590091) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(-2.5100978) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8780958) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(2.2357975) q[1];
x q[2];
rz(-1.8144242) q[3];
sx q[3];
rz(-1.676179) q[3];
sx q[3];
rz(-1.6107314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(1.9155496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292424) q[0];
sx q[0];
rz(-1.3717522) q[0];
sx q[0];
rz(-2.6890486) q[0];
rz(-pi) q[1];
rz(-2.366757) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(1.5409842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5495758) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(-1.7906584) q[1];
x q[2];
rz(-0.99818228) q[3];
sx q[3];
rz(-1.5787573) q[3];
sx q[3];
rz(1.5919459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(1.0294611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80699608) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(-2.8732804) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6497518) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(2.8284555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2468977) q[1];
sx q[1];
rz(-0.20935911) q[1];
sx q[1];
rz(0.73114242) q[1];
rz(-pi) q[2];
rz(-3.0089278) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.9784137) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.5159336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3665873) q[0];
sx q[0];
rz(-0.89582743) q[0];
sx q[0];
rz(1.3875899) q[0];
rz(-pi) q[1];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(2.5474472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.792946) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(-1.2299041) q[1];
rz(-pi) q[2];
rz(2.5528615) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(-2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-0.89486665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1631854) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(1.0582256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.52389) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(0.4862116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-1.1333229) q[1];
sx q[1];
rz(0.84769627) q[1];
x q[2];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.8572516) q[3];
sx q[3];
rz(-2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(-1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(2.6766052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156292) q[0];
sx q[0];
rz(-2.0032126) q[0];
sx q[0];
rz(0.36884357) q[0];
x q[1];
rz(-2.8000185) q[2];
sx q[2];
rz(-1.5014868) q[2];
sx q[2];
rz(2.9937033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(2.2062917) q[1];
rz(-2.898071) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(3.0628672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-2.9856317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10931817) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(-0.59662915) q[0];
rz(-pi) q[1];
rz(-1.9975125) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(2.0063426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(-3.0867982) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74219269) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.5163039) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(-0.46778932) q[0];
rz(-1.0609264) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(2.5271497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(-2.2680125) q[1];
x q[2];
rz(-1.1484654) q[3];
sx q[3];
rz(-2.9557807) q[3];
sx q[3];
rz(-1.2687792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(1.1286939) q[0];
x q[1];
rz(0.73762383) q[2];
sx q[2];
rz(-0.80543033) q[2];
sx q[2];
rz(1.1951624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5732167) q[1];
sx q[1];
rz(-1.7163367) q[1];
sx q[1];
rz(0.88370609) q[1];
rz(-1.9570458) q[3];
sx q[3];
rz(-2.2830314) q[3];
sx q[3];
rz(-1.4854747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.261196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226782) q[0];
sx q[0];
rz(-2.7180674) q[0];
sx q[0];
rz(-0.25869297) q[0];
x q[1];
rz(1.5485498) q[2];
sx q[2];
rz(-2.1689479) q[2];
sx q[2];
rz(-1.0215789) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0295769) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(-0.13346787) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19461467) q[3];
sx q[3];
rz(-0.6932887) q[3];
sx q[3];
rz(0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5579055) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-1.530004) q[2];
sx q[2];
rz(-0.57208021) q[2];
sx q[2];
rz(-0.013442599) q[2];
rz(2.9363587) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
