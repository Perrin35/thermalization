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
rz(-2.0071964) q[0];
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
rz(1.3463319) q[0];
sx q[0];
rz(-1.594461) q[0];
sx q[0];
rz(-2.6803826) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3590091) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(-2.5100978) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2156148) q[1];
sx q[1];
rz(-1.8324592) q[1];
sx q[1];
rz(-1.91933) q[1];
rz(-pi) q[2];
rz(-1.9840368) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81235028) q[0];
sx q[0];
rz(-1.3717522) q[0];
sx q[0];
rz(0.45254405) q[0];
rz(-0.77483564) q[2];
sx q[2];
rz(-0.92445395) q[2];
sx q[2];
rz(-1.6006084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1184428) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(-1.5688194) q[1];
x q[2];
rz(1.5854884) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(-0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(1.0294611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3345966) q[0];
sx q[0];
rz(-0.96141978) q[0];
sx q[0];
rz(-2.8732804) q[0];
rz(-pi) q[1];
rz(-1.2842032) q[2];
sx q[2];
rz(-2.2006052) q[2];
sx q[2];
rz(1.0857925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50476915) q[1];
sx q[1];
rz(-1.7261337) q[1];
sx q[1];
rz(1.4298646) q[1];
rz(-pi) q[2];
x q[2];
rz(1.449552) q[3];
sx q[3];
rz(-0.83449927) q[3];
sx q[3];
rz(-2.8672225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31947485) q[0];
sx q[0];
rz(-1.7135156) q[0];
sx q[0];
rz(-0.6832173) q[0];
x q[1];
rz(2.1298725) q[2];
sx q[2];
rz(-1.6435197) q[2];
sx q[2];
rz(-2.5474472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1277395) q[1];
sx q[1];
rz(-1.2989559) q[1];
sx q[1];
rz(-0.66687648) q[1];
rz(-pi) q[2];
rz(-1.3514148) q[3];
sx q[3];
rz(-2.1483768) q[3];
sx q[3];
rz(0.95336174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-0.98658371) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653523) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(2.4549237) q[0];
rz(-pi) q[1];
rz(0.61770265) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(-0.4862116) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.289031) q[1];
sx q[1];
rz(-2.3173996) q[1];
sx q[1];
rz(-2.186071) q[1];
x q[2];
rz(-0.026167913) q[3];
sx q[3];
rz(-0.28655616) q[3];
sx q[3];
rz(-2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-2.6766052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2705921) q[0];
sx q[0];
rz(-0.56068476) q[0];
sx q[0];
rz(0.90765783) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34157413) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(0.14788936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(-0.93530099) q[1];
rz(-pi) q[2];
rz(0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(2.9856317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6316846) q[0];
sx q[0];
rz(-1.317306) q[0];
sx q[0];
rz(-1.7476728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1440802) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(-1.13525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29282001) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(-1.8715026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0714674) q[3];
sx q[3];
rz(-2.0519749) q[3];
sx q[3];
rz(-1.037998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(-2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45390689) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.2841671) q[0];
rz(2.3209004) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12411815) q[1];
sx q[1];
rz(-2.0548471) q[1];
sx q[1];
rz(2.2494621) q[1];
rz(-pi) q[2];
rz(-1.1484654) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(1.2687792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(3.0548813) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(2.8588262) q[0];
rz(-2.4400318) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(-1.823002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9160999) q[0];
sx q[0];
rz(-0.44888228) q[0];
sx q[0];
rz(1.7569957) q[0];
x q[1];
rz(0.65638541) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(0.93751794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(-1.3436951) q[1];
x q[2];
rz(2.7301844) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5426853) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.8803966) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0364089) q[0];
sx q[0];
rz(-1.9793708) q[0];
sx q[0];
rz(-1.4559792) q[0];
rz(-pi) q[1];
rz(-2.5433259) q[2];
sx q[2];
rz(-1.5524128) q[2];
sx q[2];
rz(-0.53668864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65834261) q[1];
sx q[1];
rz(-1.5066506) q[1];
sx q[1];
rz(0.49883962) q[1];
rz(1.4114755) q[3];
sx q[3];
rz(-0.89309249) q[3];
sx q[3];
rz(-0.89980984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(2.5718001) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8284843) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(-2.5333511) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-0.99909487) q[2];
sx q[2];
rz(-1.5487164) q[2];
sx q[2];
rz(-1.6185417) q[2];
rz(2.5504997) q[3];
sx q[3];
rz(-1.4553087) q[3];
sx q[3];
rz(1.6607264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
