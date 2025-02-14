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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(2.5084578) q[1];
sx q[1];
rz(8.8517744) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6010701) q[0];
sx q[0];
rz(-0.94694505) q[0];
sx q[0];
rz(0.076657587) q[0];
rz(2.6296205) q[2];
sx q[2];
rz(-2.3880771) q[2];
sx q[2];
rz(-1.7640424) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4595653) q[1];
sx q[1];
rz(-0.50217705) q[1];
sx q[1];
rz(-0.14634303) q[1];
x q[2];
rz(1.6735733) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(2.4569907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28213349) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(2.0852883) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(-1.4200042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2317155) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(2.7114482) q[0];
rz(3.0138956) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970222) q[0];
sx q[0];
rz(-1.26957) q[0];
sx q[0];
rz(1.3006849) q[0];
x q[1];
rz(3.0906648) q[2];
sx q[2];
rz(-1.6279334) q[2];
sx q[2];
rz(-2.6661851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.592652) q[1];
sx q[1];
rz(-2.2363064) q[1];
sx q[1];
rz(1.1489465) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53967472) q[3];
sx q[3];
rz(-1.3840578) q[3];
sx q[3];
rz(0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0251856) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(0.44542584) q[2];
rz(-0.40924254) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55260783) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(0.71480042) q[0];
rz(2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(2.8143299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11571685) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(2.9469423) q[0];
rz(0.43250044) q[2];
sx q[2];
rz(-1.0446781) q[2];
sx q[2];
rz(-1.8217349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63673692) q[1];
sx q[1];
rz(-1.5914006) q[1];
sx q[1];
rz(-1.7223174) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9182358) q[3];
sx q[3];
rz(-1.8104324) q[3];
sx q[3];
rz(-2.6289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0032234) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.9085931) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656146) q[0];
sx q[0];
rz(-0.28802832) q[0];
sx q[0];
rz(-2.0095129) q[0];
rz(-pi) q[1];
rz(0.92278752) q[2];
sx q[2];
rz(-1.4013259) q[2];
sx q[2];
rz(-0.10709912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.33267) q[1];
sx q[1];
rz(-2.4415209) q[1];
sx q[1];
rz(-1.1247404) q[1];
rz(-1.5855012) q[3];
sx q[3];
rz(-1.8222408) q[3];
sx q[3];
rz(-2.3594643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(-2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(-2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.1530217) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(-2.8797454) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76806289) q[0];
sx q[0];
rz(-3.0515915) q[0];
sx q[0];
rz(-1.1885719) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1543617) q[2];
sx q[2];
rz(-1.2886184) q[2];
sx q[2];
rz(-0.63167494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3186372) q[1];
sx q[1];
rz(-1.3483817) q[1];
sx q[1];
rz(0.59808029) q[1];
rz(-2.35942) q[3];
sx q[3];
rz(-1.7165136) q[3];
sx q[3];
rz(0.71469939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(0.74553472) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(-0.44031269) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(0.17131677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765613) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(1.2622467) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8181658) q[2];
sx q[2];
rz(-0.77980552) q[2];
sx q[2];
rz(-0.16743539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3503987) q[1];
sx q[1];
rz(-1.5916675) q[1];
sx q[1];
rz(-0.69857614) q[1];
rz(0.8796575) q[3];
sx q[3];
rz(-0.98239693) q[3];
sx q[3];
rz(-0.0042875687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(0.94179955) q[2];
rz(-1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.340516) q[3];
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
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-2.4536224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82179994) q[0];
sx q[0];
rz(-1.742779) q[0];
sx q[0];
rz(-0.71806192) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9161467) q[2];
sx q[2];
rz(-0.98042578) q[2];
sx q[2];
rz(1.2638059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35536218) q[1];
sx q[1];
rz(-1.3578102) q[1];
sx q[1];
rz(-2.0128241) q[1];
rz(2.4554208) q[3];
sx q[3];
rz(-2.1451575) q[3];
sx q[3];
rz(-2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(-0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(-3.0514858) q[0];
rz(1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(2.7065014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6259436) q[0];
sx q[0];
rz(-1.0727779) q[0];
sx q[0];
rz(1.7367212) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7516364) q[2];
sx q[2];
rz(-2.2308084) q[2];
sx q[2];
rz(1.7969799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96174091) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(-2.1895616) q[1];
rz(-pi) q[2];
rz(0.42436142) q[3];
sx q[3];
rz(-1.5699982) q[3];
sx q[3];
rz(0.80959807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(2.5153861) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597252) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(-3.0873121) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.3351701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568118) q[0];
sx q[0];
rz(-0.81997141) q[0];
sx q[0];
rz(0.7818082) q[0];
x q[1];
rz(2.1303875) q[2];
sx q[2];
rz(-1.8768684) q[2];
sx q[2];
rz(2.0843506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74354913) q[1];
sx q[1];
rz(-0.7760007) q[1];
sx q[1];
rz(0.77568027) q[1];
rz(-2.0306251) q[3];
sx q[3];
rz(-2.9879192) q[3];
sx q[3];
rz(0.75017649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(0.35150251) q[2];
rz(1.3737804) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(-0.75916284) q[0];
rz(-2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-1.0677451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7368374) q[0];
sx q[0];
rz(-0.99946293) q[0];
sx q[0];
rz(-0.57147567) q[0];
rz(1.8482089) q[2];
sx q[2];
rz(-1.2870064) q[2];
sx q[2];
rz(-2.0744155) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2764923) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(-0.92030763) q[1];
rz(-pi) q[2];
rz(2.8078733) q[3];
sx q[3];
rz(-0.90147831) q[3];
sx q[3];
rz(-3.1143318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4296253) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(2.9313226) q[2];
rz(-2.353904) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.632985) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(-1.8456712) q[3];
sx q[3];
rz(-1.4209005) q[3];
sx q[3];
rz(3.0358547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
