OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4700692) q[0];
sx q[0];
rz(-2.1052261) q[0];
sx q[0];
rz(-0.1202017) q[0];
rz(-1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(-2.3766975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51812664) q[1];
sx q[1];
rz(-1.4033485) q[1];
sx q[1];
rz(-1.8159588) q[1];
x q[2];
rz(-2.8005881) q[3];
sx q[3];
rz(-1.4031646) q[3];
sx q[3];
rz(2.117702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6630163) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-0.93634161) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124356) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-0.23075128) q[0];
x q[1];
rz(0.41994862) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(0.74479693) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(-2.9794934) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6132658) q[3];
sx q[3];
rz(-0.50308933) q[3];
sx q[3];
rz(-0.79723061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24580978) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9868601) q[0];
sx q[0];
rz(-1.2504471) q[0];
sx q[0];
rz(-2.8264168) q[0];
rz(-pi) q[1];
rz(2.2592696) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(-2.4633212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(1.1413241) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9150312) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(-1.759474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216953) q[0];
sx q[0];
rz(-1.2384402) q[0];
sx q[0];
rz(-0.38145782) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5449764) q[2];
sx q[2];
rz(-0.46586793) q[2];
sx q[2];
rz(-2.4398838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(-2.1898502) q[1];
rz(-pi) q[2];
rz(-1.0700978) q[3];
sx q[3];
rz(-0.90161937) q[3];
sx q[3];
rz(-2.5239528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-1.0505189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8331063) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(-1.2480877) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6987223) q[2];
sx q[2];
rz(-1.0204698) q[2];
sx q[2];
rz(1.0843104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1615636) q[1];
sx q[1];
rz(-1.0161576) q[1];
sx q[1];
rz(2.5142923) q[1];
rz(-pi) q[2];
rz(-2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(-2.7725078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(-1.9472306) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69960064) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(1.4621428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7397241) q[0];
sx q[0];
rz(-1.3603633) q[0];
sx q[0];
rz(1.4896638) q[0];
x q[1];
rz(0.92163779) q[2];
sx q[2];
rz(-1.7643133) q[2];
sx q[2];
rz(-2.9606539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6979132) q[1];
sx q[1];
rz(-1.7410478) q[1];
sx q[1];
rz(2.7531569) q[1];
rz(-pi) q[2];
rz(-1.9147647) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(-2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2016466) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-2.3366826) q[2];
rz(1.9645875) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(-2.129508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18400684) q[0];
sx q[0];
rz(-2.036096) q[0];
sx q[0];
rz(2.1561949) q[0];
x q[1];
rz(0.38883932) q[2];
sx q[2];
rz(-1.6992339) q[2];
sx q[2];
rz(1.5156137) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.583657) q[1];
sx q[1];
rz(-2.5968938) q[1];
sx q[1];
rz(-0.063636585) q[1];
x q[2];
rz(3.0104962) q[3];
sx q[3];
rz(-0.56250611) q[3];
sx q[3];
rz(-1.2870115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(-1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-0.28373757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34127125) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(-0.11462258) q[0];
rz(-1.2405422) q[2];
sx q[2];
rz(-1.7527765) q[2];
sx q[2];
rz(0.44588003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(-2.6372361) q[1];
rz(-0.65225668) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(3.101845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5725296) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208961) q[0];
sx q[0];
rz(-0.81917742) q[0];
sx q[0];
rz(0.92185123) q[0];
rz(2.9022129) q[2];
sx q[2];
rz(-2.8068672) q[2];
sx q[2];
rz(-0.3604381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3410586) q[1];
sx q[1];
rz(-1.3375999) q[1];
sx q[1];
rz(-1.484616) q[1];
rz(-pi) q[2];
rz(1.2980372) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(-0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(-0.64176732) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(0.26783255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(-2.3964336) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7540625) q[2];
sx q[2];
rz(-1.6943309) q[2];
sx q[2];
rz(-1.087041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99991998) q[1];
sx q[1];
rz(-1.281764) q[1];
sx q[1];
rz(2.8029867) q[1];
x q[2];
rz(-2.0249428) q[3];
sx q[3];
rz(-2.5513259) q[3];
sx q[3];
rz(2.2772307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(-2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.98942479) q[2];
sx q[2];
rz(-0.94716723) q[2];
sx q[2];
rz(-0.92826044) q[2];
rz(0.6775425) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
